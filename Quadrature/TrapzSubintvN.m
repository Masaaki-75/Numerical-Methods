function val = TrapzSubintvN(x,f,nintv)
%==========================================================================
%TRAPZSUBINTVN return the numerical quadrature using Trapezoidal method.
%
% Author: Ma Chenglong (github: Masaaki-75)
%         Copyright 2021 Sun Yat-sen University.
%
% Syntax:
%
%   VAL = TRAPZSUBINTVN(X,F,NINTV);
%       
%       X is a vector containing points in the interval of the integral; F
%       is a function handle for integrand; NINTV is the number of
%       subintervals.
%
%==========================================================================

f = convertStringsToChars(f);
if ~exist('nintv','var')||isempty(nintv), nintv = 1; end
if isa(f,'sym')||isa(f,'symfun'), f = matlabFunction(f); end
if isa(f,'char'), f = eval(['@(x)',f]); end
%if isa(y,'function_handle'), y = y(x); end

if nintv > 1
    % insert (NINTV-1) points between every 2 points in x.
    xi = NInterpBtwn2(x,nintv-1);
else
    xi = x;
end

if isnumeric(f)
    yi = f;
else
    yi = f(xi);
end

if ~all(size(xi)==size(yi))
    xi = convertToVec(xi,'row');
    yi = convertToVec(yi,'row');
    if ~all(size(xi)==size(yi))
        error('Error: Input X and Y have different sizes.')
    end
end

val = trapz(xi,yi);

end

