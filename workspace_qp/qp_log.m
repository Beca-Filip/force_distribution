classdef qp_log < handle
%QP_LOG  Crash-proof run log for the long IO fits.
%
%   lg = QP_LOG(results_dir)
%   lg = QP_LOG(results_dir, 'Append', true)
%
%   A fit takes the better part of a day on a remote workstation, and the
%   thing that goes wrong is never the thing you were watching.  The second
%   B5 attempt was lost twice over: once when a CasADi error killed the batch
%   at the second condition, and again when the remote desktop rebooted and
%   took the scrollback with it.  Everything this class does follows from
%   that: write to disk, write immediately, and write enough that the run can
%   be diagnosed after the fact by someone who was not watching.
%
%
%   WHAT IT WRITES, all under results_dir
%   -------------------------------------
%     log.txt          the human stream.  Every line timestamped, every
%                      section banner-delimited.  Also carries the fmincon
%                      'iter' table, via the diary.
%     events.jsonl     one JSON object per line, one line per event.  This is
%                      the machine-readable half: what SUMMARIZE_MAIN_RESULTS
%                      and any later post-processing should read, because it
%                      does not need parsing rules that change whenever a
%                      printf is edited.
%     iterations.csv   the fmincon trace -- condition, iteration, fval,
%                      first-order optimality, step size, constraint
%                      violation, elapsed seconds.  Written per iteration, so
%                      the shape of a stalled run is visible while it is
%                      still running.
%     status.txt       ONE line, overwritten: where the run is right now and
%                      how long it has been there.  Small enough to read over
%                      a laggy remote session without opening anything.
%     checkpoint.mat   the latest theta, refreshed every fmincon iteration.
%
%   Every write opens, writes and closes the file.  MATLAB gives no way to
%   flush a held file handle, and a handle held across a power cut loses
%   whatever was buffered -- which is precisely the case this class exists
%   for.  The cost is a few file operations per iteration against multi-second
%   iterations, so it does not register.
%
%
%   TYPICAL USE
%   -----------
%       lg = qp_log(results_dir);
%       lg.header('B5 per-condition batch', cfg);
%       lg.section('Subject 4  speed 1  leg 2');
%       lg.status('subject 4 speed 1 leg 2: fitting');
%       lg.printf('cond(Q) = %.3e', cond_Q);
%       lg.event('qp_degraded', struct('sample', k, 'rung', info.rung));
%       ...
%       lg.exception(ME, struct('subject', 4, 'speed', 1, 'leg', 2));
%       lg.close();
%
%   The fmincon hook is
%       'OutputFcn', @(th, ov, st) lg.fmincon_hook(th, ov, st, tag)
%   which writes an iterations.csv row and refreshes checkpoint.mat.  It
%   never throws: a logging failure must not end a fit that is going fine.
%
%   See also MAIN, QP_HEALTH_CHECK, QP_SOLVE, SUMMARIZE_MAIN_RESULTS.

    properties (SetAccess = private)
        dir             % results directory
        log_path        % log.txt
        events_path     % events.jsonl
        iters_path      % iterations.csv
        status_path     % status.txt
        ckpt_path       % checkpoint.mat
        t_start         % tic at construction
        n_warn = 0      % counters, reported by close()
        n_error = 0
    end

    methods
        function obj = qp_log(results_dir, varargin)
            p = inputParser;
            p.addParameter('Append', false, @(x) islogical(x) || isnumeric(x));
            p.addParameter('Diary',  true,  @(x) islogical(x) || isnumeric(x));
            p.parse(varargin{:});

            if ~exist(results_dir, 'dir')
                mkdir(results_dir);
            end

            obj.dir         = results_dir;
            obj.log_path    = fullfile(results_dir, 'log.txt');
            obj.events_path = fullfile(results_dir, 'events.jsonl');
            obj.iters_path  = fullfile(results_dir, 'iterations.csv');
            obj.status_path = fullfile(results_dir, 'status.txt');
            obj.ckpt_path   = fullfile(results_dir, 'checkpoint.mat');
            obj.t_start     = tic;

            if ~p.Results.Append
                % Roll the previous run aside rather than overwrite it.  A
                % rerun almost always follows a failure, and the failed run's
                % log is the thing you want to compare against.
                stamp = datestr(now, 'yyyymmdd-HHMMSS');  %#ok<TNOW1,DATST>
                for f = {obj.log_path, obj.events_path, obj.iters_path}
                    if exist(f{1}, 'file')
                        [d, b, e] = fileparts(f{1});
                        movefile(f{1}, fullfile(d, sprintf('%s-%s%s', b, stamp, e)));
                    end
                end
            end

            if ~exist(obj.iters_path, 'file')
                obj.append(obj.iters_path, sprintf( ...
                    'condition,iteration,funccount,fval,firstorderopt,stepsize,constrviolation,elapsed_s\n'));
            end

            if p.Results.Diary
                diary(obj.log_path);
                diary on;
            end
        end

        % -----------------------------------------------------------------
        function header(obj, title, cfg)
            %HEADER  Banner plus the provenance a post-mortem needs.
            %
            %   Which code produced a result matters more here than usual:
            %   the first B5 batch was audited months later and turned out to
            %   have been run from a commit that predated the conditioning
            %   work, which was only provable from what the .mat files did
            %   NOT contain.  Recording the commit up front makes that a
            %   lookup instead of an inference.
            obj.rule('=');
            obj.printf_raw(' %s', title);
            obj.printf_raw(' started      %s', datestr(now, 'yyyy-mm-dd HH:MM:SS'));  %#ok<TNOW1,DATST>
            obj.printf_raw(' host         %s', qp_log.hostname());
            obj.printf_raw(' matlab       %s', version());
            obj.printf_raw(' git commit   %s', qp_log.git_commit());
            obj.printf_raw(' output       %s', obj.dir);
            obj.rule('=');

            if nargin > 2 && ~isempty(cfg)
                obj.printf_raw('Configuration:');
                obj.print_struct(cfg, '  ');
                obj.rule('-');
            end

            obj.event('run_start', struct( ...
                'title',   title, ...
                'host',    qp_log.hostname(), ...
                'matlab',  version(), ...
                'commit',  qp_log.git_commit()));
        end

        function section(obj, fmt, varargin)
            obj.rule('-');
            obj.printf_raw([' ' fmt], varargin{:});
            obj.rule('-');
        end

        function printf(obj, fmt, varargin)
            %PRINTF  A timestamped line to the console and to log.txt.
            line = sprintf(fmt, varargin{:});
            fprintf('[%s] %s\n', obj.clock(), line);
            obj.mirror(sprintf('[%s] %s\n', obj.clock(), line));
        end

        function printf_raw(obj, fmt, varargin)
            %PRINTF_RAW  As printf, without the timestamp.  For banners.
            line = sprintf(fmt, varargin{:});
            fprintf('%s\n', line);
            obj.mirror(sprintf('%s\n', line));
        end

        function warn(obj, id, fmt, varargin)
            obj.n_warn = obj.n_warn + 1;
            msg = sprintf(fmt, varargin{:});
            obj.printf('WARNING (%s) %s', id, msg);
            obj.event('warning', struct('id', id, 'message', msg));
        end

        % -----------------------------------------------------------------
        function status(obj, fmt, varargin)
            %STATUS  Overwrite the one-line "where is it now" file.
            msg = sprintf(fmt, varargin{:});
            txt = sprintf('%s | elapsed %s | %s\n', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'), ...
                qp_log.hms(toc(obj.t_start)), msg);  %#ok<TNOW1,DATST>
            obj.overwrite(obj.status_path, txt);
        end

        function event(obj, kind, payload)
            %EVENT  Append one JSON object to events.jsonl.
            %
            %   Logging must never be the reason a fit dies, so a payload
            %   that will not encode is reported as a string rather than
            %   raised.
            if nargin < 3
                payload = struct();
            end
            rec = struct( ...
                't',         datestr(now, 'yyyy-mm-ddTHH:MM:SS'), ...  %#ok<TNOW1,DATST>
                'elapsed_s', toc(obj.t_start), ...
                'kind',      kind, ...
                'data',      payload);
            try
                txt = jsonencode(rec);
            catch
                txt = jsonencode(struct('t', rec.t, 'elapsed_s', rec.elapsed_s, ...
                    'kind', kind, 'data', '<payload not JSON-encodable>'));
            end
            obj.append(obj.events_path, [txt sprintf('\n')]);
        end

        function exception(obj, ME, context)
            %EXCEPTION  The full story of a failure, in both files.
            %
            %   Message AND identifier AND stack.  The identifier is what
            %   later code branches on, the stack is what a reader needs, and
            %   CasADi's message body is a screenful of raw solver inputs
            %   that is worth keeping but not worth reading first -- so it
            %   goes last.
            if nargin < 3
                context = struct();
            end
            obj.n_error = obj.n_error + 1;

            obj.rule('!');
            obj.printf_raw(' FAILED: %s', ME.identifier);
            if ~isempty(fieldnames(context))
                obj.printf_raw(' context:');
                obj.print_struct(context, '   ');
            end
            obj.printf_raw(' stack:');
            for s = 1:numel(ME.stack)
                obj.printf_raw('   %s (line %d)', ME.stack(s).name, ME.stack(s).line);
            end
            obj.printf_raw(' message:');
            msg_lines = strsplit(ME.message, newline);
            for s = 1:min(numel(msg_lines), 20)
                obj.printf_raw('   %s', msg_lines{s});
            end
            if numel(msg_lines) > 20
                obj.printf_raw('   ... %d more lines, full text in events.jsonl', ...
                    numel(msg_lines) - 20);
            end
            obj.rule('!');

            stack = struct('name', {}, 'line', {});
            for s = 1:numel(ME.stack)
                stack(s).name = ME.stack(s).name;
                stack(s).line = ME.stack(s).line;
            end
            obj.event('exception', struct( ...
                'identifier', ME.identifier, ...
                'message',    ME.message, ...
                'stack',      stack, ...
                'context',    context));
        end

        % -----------------------------------------------------------------
        function stop = fmincon_hook(obj, theta, optimValues, state, tag)
            %FMINCON_HOOK  OutputFcn: one CSV row and one checkpoint per iteration.
            %
            %   Wrapped whole in try/catch.  An OutputFcn that throws aborts
            %   fmincon, and a run must not be lost because a disk was busy.
            stop = false;
            if nargin < 5
                tag = 'fit';
            end
            try
                if strcmp(state, 'iter')
                    obj.append(obj.iters_path, sprintf('%s,%d,%d,%.10g,%.10g,%.10g,%.10g,%.3f\n', ...
                        tag, ...
                        optimValues.iteration, ...
                        qp_log.field_or(optimValues, 'funccount', NaN), ...
                        optimValues.fval, ...
                        qp_log.field_or(optimValues, 'firstorderopt', NaN), ...
                        qp_log.field_or(optimValues, 'stepsize', NaN), ...
                        qp_log.field_or(optimValues, 'constrviolation', NaN), ...
                        toc(obj.t_start)));

                    ckpt = struct( ...
                        'tag',             tag, ...
                        'theta',           theta, ...
                        'iteration',       optimValues.iteration, ...
                        'fval',            optimValues.fval, ...
                        'firstorderopt',   qp_log.field_or(optimValues, 'firstorderopt', NaN), ...
                        'constrviolation', qp_log.field_or(optimValues, 'constrviolation', NaN), ...
                        'elapsed_s',       toc(obj.t_start), ...
                        'saved_at',        datestr(now, 'yyyy-mm-dd HH:MM:SS'));  %#ok<TNOW1,DATST>
                    save(obj.ckpt_path, '-struct', 'ckpt');

                    obj.status('%s: iteration %d, fval %.6g, firstorderopt %.3g', ...
                        tag, optimValues.iteration, optimValues.fval, ...
                        qp_log.field_or(optimValues, 'firstorderopt', NaN));
                end
            catch
                % Deliberately silent: reporting a logging failure through
                % the log is circular, and the fit is what matters.
            end
        end

        % -----------------------------------------------------------------
        function close(obj)
            obj.rule('=');
            obj.printf_raw(' finished %s after %s   (%d warnings, %d failures)', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'), ...
                qp_log.hms(toc(obj.t_start)), obj.n_warn, obj.n_error);  %#ok<TNOW1,DATST>
            obj.rule('=');
            obj.event('run_end', struct('elapsed_s', toc(obj.t_start), ...
                'n_warn', obj.n_warn, 'n_error', obj.n_error));
            obj.status('finished: %d warnings, %d failures', obj.n_warn, obj.n_error);
            diary('off');
        end
    end

    % =====================================================================
    methods (Access = private)
        function rule(obj, ch)
            obj.printf_raw('%s', repmat(ch, 1, 73));
        end

        function print_struct(obj, s, indent)
            names = fieldnames(s);
            for i = 1:numel(names)
                obj.printf_raw('%s%-16s %s', indent, names{i}, ...
                    qp_log.brief(s.(names{i})));
            end
        end

        function mirror(obj, txt)
            % The diary already mirrors the console into log.txt.  Without a
            % diary (a worker, a nested call) nothing would reach the file,
            % so write it explicitly in that case only, to avoid doubling
            % every line.
            if strcmp(get(0, 'Diary'), 'off')
                obj.append(obj.log_path, txt);
            end
        end

        function append(~, path, txt)
            fid = fopen(path, 'a');
            if fid > 0
                fprintf(fid, '%s', txt);
                fclose(fid);
            end
        end

        function overwrite(~, path, txt)
            fid = fopen(path, 'w');
            if fid > 0
                fprintf(fid, '%s', txt);
                fclose(fid);
            end
        end

        function s = clock(~)
            s = datestr(now, 'HH:MM:SS');  %#ok<TNOW1,DATST>
        end
    end

    % =====================================================================
    methods (Static, Access = private)
        function v = field_or(s, name, default)
            if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
                v = s.(name);
            else
                v = default;
            end
        end

        function h = hostname()
            h = getenv('COMPUTERNAME');
            if isempty(h)
                h = getenv('HOSTNAME');
            end
            if isempty(h)
                h = 'unknown';
            end
        end

        function c = git_commit()
            %GIT_COMMIT  Short hash plus a dirty marker, or 'unavailable'.
            c = 'unavailable';
            here = fileparts(mfilename('fullpath'));
            try
                [st, out] = system(sprintf('git -C "%s" rev-parse --short HEAD', here));
                if st == 0
                    c = strtrim(out);
                    [st2, out2] = system(sprintf('git -C "%s" status --porcelain', here));
                    if st2 == 0 && ~isempty(strtrim(out2))
                        c = [c ' (dirty)'];
                    end
                end
            catch
                % leave 'unavailable'
            end
        end

        function s = hms(t)
            s = sprintf('%02d:%02d:%02d', floor(t/3600), ...
                mod(floor(t/60), 60), mod(floor(t), 60));
        end

        function s = brief(v)
            %BRIEF  One-line rendering of a config value.
            if ischar(v) || isstring(v)
                s = char(v);
            elseif islogical(v)
                s = mat2str(v);
            elseif isnumeric(v) && isscalar(v)
                s = sprintf('%g', v);
            elseif isnumeric(v) && numel(v) <= 12
                s = mat2str(v(:)', 6);
            elseif isnumeric(v)
                s = sprintf('[%s %s], %g .. %g', mat2str(size(v)), class(v), ...
                    min(v(:)), max(v(:)));
            else
                s = sprintf('<%s %s>', mat2str(size(v)), class(v));
            end
        end
    end
end
