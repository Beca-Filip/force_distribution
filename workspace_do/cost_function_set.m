function J = cost_function_set(f, fmin, fmax, pcsa, vmt, M, fpassive, f0, m, r)

% J = zeros(1, 17, 'like', casadi.MX);
% 
% J(1) = cost_function1(f);
% J(2) = cost_function2(f);
% J(3) = cost_function3(f);
% J(4) = cost_function4(f);
% J(5) = cost_function5(f, fmin, fmax);
% J(6) = cost_function6(f, fmin, fmax);
% J(7) = cost_function7(f, fmin, fmax);
% J(8) = cost_function8(f, fmin, fmax);
% J(9) = cost_function9(f, pcsa);
% J(10) = cost_function10(f, pcsa);
% J(11) = cost_function11(f, pcsa);
% J(12) = cost_function12(f, pcsa);
% J(13) = cost_function13(f, vmt);
% J(14) = cost_function14(f, M);
% J(15) = cost_function15(f, fpassive, f0, pcsa, m);
% J(16) = cost_function16(f, f0, r);
% J(17) = cost_function17(f, fmax);
% 
% excl = [2:17];
% Jnew = zeros(1, 17 - length(excl), 'like', casadi.MX);
% cnt = 1;
% for ii = 1 : 17
%     if ~ismember(ii, excl)
%         Jnew(cnt) = J(ii);
%         cnt = cnt+1;
%     end
% end
% J = Jnew;

excl = [17];
J = zeros(1, 17 - length(excl), 'like', casadi.MX);


cnt = 1;
if ~ismember(1, excl)
J(cnt) = cost_function1(f);
cnt = cnt+1;
end

if ~ismember(2, excl)
J(cnt) = cost_function2(f);
cnt = cnt+1;
end

if ~ismember(3, excl)
J(cnt) = cost_function3(f);
cnt = cnt+1;
end

if ~ismember(4, excl)
J(cnt) = cost_function4(f);
cnt = cnt+1;
end

if ~ismember(5, excl)
J(cnt) = cost_function5(f, fmin, fmax);
cnt = cnt+1;
end

if ~ismember(6, excl)
J(cnt) = cost_function6(f, fmin, fmax);
cnt = cnt+1;
end

if ~ismember(7, excl)
J(cnt) = cost_function7(f, fmin, fmax);
cnt = cnt+1;
end

if ~ismember(8, excl)
J(cnt) = cost_function8(f, fmin, fmax);
cnt = cnt+1;
end

if ~ismember(9, excl)
J(cnt) = cost_function9(f, pcsa);
cnt = cnt+1;
end

if ~ismember(10, excl)
J(cnt) = cost_function10(f, pcsa);
cnt = cnt+1;
end

if ~ismember(11, excl)
J(cnt) = cost_function11(f, pcsa);
cnt = cnt+1;
end

if ~ismember(12, excl)
J(cnt) = cost_function12(f, pcsa);
cnt = cnt+1;
end

if ~ismember(13, excl)
J(cnt) = cost_function13(f, vmt);
cnt = cnt+1;
end

if ~ismember(14, excl)
J(cnt) = cost_function14(f, M);
cnt = cnt+1;
end

if ~ismember(15, excl)
J(cnt) = cost_function15(f, fpassive, f0, pcsa, m);
cnt = cnt+1;
end

if ~ismember(16, excl)
J(cnt) = cost_function16(f, f0, r);
cnt = cnt+1;
end

if ~ismember(17, excl)
J(cnt) = cost_function17(f, fmax);
cnt = cnt+1;
end


end