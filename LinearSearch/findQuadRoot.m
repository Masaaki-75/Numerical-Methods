function varargout = findQuadRoot(varargin)
%==========================================================================
%FINDQUADROOT return the roots of given quadratic equation.
%
% Author: Ma Chenglong (github: Masaaki-75)
%         Copyright 2021 Sun Yat-sen University.
%
% Syntax:
%
%   X = FINDQUADROOT(A,B,C);
%
%       returns the estimated root X of the quadratic equation specified by
%       the coefficients A,B,C.
%       
%       If A,B,C are all scalars, then the equation would be
%
%           A*x^2 + B*x + C = 0;
%
%       But if A,B,C are all (column) vectors with the same length N, then
%       each element of them forms a triplet
%
%           [A(i), B(i), C(i)],   i=1,2,...,N
%
%       that depicts a quadratic equation. Therefore, the result X will be
%       a N-by-2 matrix, each row of which contains dual solution of one
%       equation. For example, when N = 3:
%
%               [a1]       [b1]      [c1]         [x11,x12]
%           A = [a2],  B = [b2], C = [c2] -> X =  [x21,x22]
%               [a3]       [b3]      [c3]         [x31,x32]
%
%
%   X = FINDQUADROOT(B,C);
%       
%       returns the estimated root X of the quadratic equation specified by
%       the coefficients B,C as follows:
%
%           x^2 + B*x + C = 0;
%
%       As mentioned above, B and C can be either scalars or vectors.
%
%   X = FINDQUADROOT(P);
%       
%       returns the estimated root X of the quadratic equation specified by
%       the coefficients A,B,C. Here, A, B and C are three columns of the
%       matrix P, i.e. 
%
%               [a1, b1, c1]          [x11,x12]
%           P = [a2, b2, c2]  ->  X = [x21,x22]
%               [a3, b3, c3]          [x31,x32]
%
%       Of course, P can also be a 3-dimension vector.
%
%
%   [X1,X2] = FINDQUADROOT(...);
%
%       returns two solutions repectively. This is equivalent to the above
%       codes:
%
%       TEMP = FINDQUADROOT(...);
%       X1 = TEMP(:,1); X2 = TEMP(:,2);
%
%==========================================================================

if nargin == 1
    p = varargin{1};
    if size(p,2) ~= 3, p = p.'; end
    validateattributes(p,{'numeric'},{'size',[NaN,3]})
    a = p(:,1); b = p(:,2); c = p(:,3);
    % [a1, b1, c1]      [x11,x12]
    % [a2, b2, c2]  ->  [x21,x22]
    % [a3, b3, c3]      [x31,x32]
elseif nargin == 2
    b = varargin{1};
    c = varargin{2};
    a = ones(size(b));
elseif nargin == 3
    a = varargin{1};
    b = varargin{2};
    c = varargin{3};
else
    error('Error: Too many input arguments.')
end

% a,b,c are all column vectors.
if size(a,2) ~= 1, a = a.'; end
if size(b,2) ~= 1, b = b.'; end
if size(c,2) ~= 1, c = c.'; end

delta = b.^2 - 4*a.*c;
sqrtd = sqrt(delta);

n_eqn = size(a,1);
x = zeros(n_eqn,2);

x(:,1) =  (-b - sqrtd)./(2*a);
x(:,2) =  (-b + sqrtd)./(2*a);
% disp(array2table([(1:size(x,1))',x],'VariableNames',{'No','x1','x2'})) 


if nargout == 2
    varargout{1} = x(:,1);
    varargout{2} = x(:,2);
else
    varargout{1} = x;
end