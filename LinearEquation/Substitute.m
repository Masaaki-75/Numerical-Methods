function x = Substitute(A,b,direction)
%==========================================================================
%SUBSTITUTE Forward or backward substitution of the linear equation system.
%
% Note that the input coefficient matrix A should be of either upper
% triangular form or lower triangular one.
%
% Author: Ma Chenglong (github: Masaaki-75)
%         Copyright 2021 Sun Yat-sen University.
%
% Syntax: 
%
%   X = SUBSTITUTE(A,B,DIR);
%
%       returns the solution X of the input linear system (depicted by
%       coefficient matrix A and constant vector B) as a column vector 
%       using substitution method. 
%
%       B should be a column vector and a row vector B will be transposed.
%
%       DIR specifies the direction of substitution. If A is an upper
%       triangular matrix, then backward substitution should be applied and
%       the value of DIR should be assigned as 'upper', 'up' or 'u'. If A
%       is lower triangular, then DIR should be 'lower', 'low' or 'l' to
%       perform forward substitution to get the right solution.
%
%==========================================================================

%% Input Arguments Processing
if size(b,1) < size(b,2), b = b.'; end
rows = length(b); cols = size(A,2);
if rows ~= cols, error('Error: The coefficient matrix should be a square matrix.'); end
if ~exist('direction','var')||isempty(direction), direction = 'upper'; end
direction = lower(convertStringsToChars(direction));

% set x as a row vector to make use of vectorized computation.
x = zeros(1,rows);

%% Substitution
switch direction
    case {'upper','up','u'}  % backward substitution
        ii_order = rows-1:-1:1;
        x(rows) = b(rows)/A(rows,end);  % initial value in the last row
        
        for ii = ii_order
            temp = sum(A(ii,ii+1:end).*x(ii+1:end));
            x(ii) =  (b(ii) - temp)/A(ii,ii);
            
            if abs(x(ii)) < 1e-10, x(ii) = 0; end
        end
    case {'lower','low','l'}  % forward substitution
        ii_order = 2:rows;
        x(1) = b(1)/A(1,1);  % initial value in the first row
        
        for ii = ii_order
            temp = sum(A(ii,1:ii-1).*x(1:ii-1));
            x(ii) = (b(ii) - temp)/A(ii,ii);
            
            if abs(x(ii)) < 1e-10, x(ii) = 0; end
        end
    otherwise
        error('Error: Argument "direction" should be either "up" or "low".')
end

x = x.';  % return x as a column vector.
end