function nA = multidim_column_norm(A, varargin)
%MULTIDIM_COLUMN_NORM calculates the vector of column norms of a matrix, or
%of a multidimensional stack of matrices.
%
%   nA = MULTIDIM_COLUMN_NORM(A)
%   nA = MULTIDIM_COLUMN_NORM(A, n)
%
%   Inputs:
%   A ~ multidimensional stack of matrices with size [n1, n2, ..., nd].
%   n ~ which norm to calculate (any positive integer or (+-) inf).
%   
%   Outputs:
%   nA ~ multidimensional stack of vectors with size [1, n2, ..., nd].

% If additional argument passed
if nargin > 1
    n = varargin{1};
else
    n = 2;
end

% If n is not an integer greater or equal to one or isn't +- inf.
if ~( n == inf || n == -inf || (n >= 1 && n == round(n)) )
   error('n should be a positive integer >= 1, inf or -inf.');
end

% Sizes of matrix stack
szA = size(A);
% Calculate the dimensionality of matrix stack
dA = length(szA);
% The size of the output vector stack (replace first dimension size with 1)
sznA = [1, szA(2:end)];
% Prealocate output
nA = nan(sznA);

% If there is only 2 dimension, return normal column norm
if dA == 2
    nA = column_norm(A, n);
% otherwise loop through dimensions
else
    
    % Loop through all other dimensions
    % Loop code at 
    % https://fr.mathworks.com/matlabcentral/answers/386638-how-to-loop-through-an-unknown-number-of-matrix-dimensions
    nv   = dA - 2;
    v    = ones(1, nv);
    vLim = szA(3:dA);
    ready = false;
    while ~ready
        %%%%% BEGIN LOOP
        % Transform index vector to cell
        cellv = num2cell(v, nv);
        % Calculate
        nA(1, :, cellv{:}) = column_norm(squeeze(A(:, :, cellv{:})), n);
        %%%%% END LOOP
        % Update the index vector:
        ready = true;       % Assume that the WHILE loop is ready
        for k = 1:nv
            v(k) = v(k) + 1;
            if v(k) <= vLim(k)
                ready = false;  % No, WHILE loop is not ready now
                break;          % v(k) increased successfully, leave "for k" loop
            end
            v(k) = 1;         % Reset v(k), proceed to next k
        end
    end
    
end