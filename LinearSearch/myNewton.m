function [x,results] = myNewton(f,x0,es,iter_max,silent)
%==========================================================================
%MYNEWTON estimate the root of univariate equation f(x)=0 using
% Newton-Raphson method.
%
% Author: Ma Chenglong (github: Masaaki-75)
%         Copyright 2021 Sun Yat-sen University.
%
% Syntax:
%
%   [X,RESULTS] = MYNEWTON(F,X0[,ES,ITER_MAX,SILENT]);
%
%       returns the estimated zero X and a matrix RESULTS storing the
%       relevant parameters computed in each iteration, given input 
%       function F, initial point X0, error tolerance ES and maximum number
%       of iteration ITER_MAX. The last argument SILENT determines whether 
%       to display the iteration results in command line window.
%
%       F is the given function, whose class can be one of the following
%       {'string','char','function_handle','sym','symfun'}. This argument
%       will be treated as a symbolic function in MYNEWTON() function.
%
%       X0 is a numeric scalar for initial point.
%
%       ES is the preset absolute relative error as error tolerance.
%       Optional, 1e-5 as default.
%
%       ITER_MAX is a positive integer representing the maximum number of
%       iteration. Optional, 200 as default.
%
%       SILENT is a logical value determining the display of the iteration
%       results. Set SILENT as false will print in the command line window
%       the results in each iteration. Optional, true as default. 
%
%==========================================================================

%% Input Arguments Processing
f = convertStringsToChars(f);
if isa(f,'function_handle'), f = sym(f); end
if isa(f,'char'), f = str2sym(f); end
if isa(f,'sym'), f = symfun(f,symvar(f)); end
% If the derivative of the function turns out to be a constant, then the 
% function handle of the derivative could not receive any input.
% For example, if the input function handle is
%     F = @(x) (x+1);
% then the function handle of the derivative is given by
%     FX = matlabFunction(diff(sym(F)));
% which yields
%     FX == @() 1.0;
% That is, FX returns 1.0 for any input number. However, in MATLAB,
% an error will arise if we enter any code like this:
%     FX(1)
% The error indicates that "There are too many input arguments."
% To avoid problems caused by differentiation on function handle, we
% convert it into a symbolic variable.

if ~exist('es','var')||isempty(es), es = 1e-5; end
if ~exist('iter_max','var')||isempty(iter_max), iter_max = 200; end
if ~exist('silent','var')||isempty(silent), silent = false; end
silent = lower(convertStringsToChars(silent));
if ~islogical(silent)
    switch silent
        case {'on','disp','display','show','false'}
            silent = false;
        otherwise
            silent = true;
    end
end
f1 = diff(f); x00 = x0;

%% Initial Setup
iter = 0; ea = 1;
f0 = f(x0);
results = zeros(iter_max+1,4);
results(1,:) = [iter,x0,f0,ea];

%% Iteration
while (ea > es) && (iter < iter_max)
    x = x0 - f0/f1(x0); 
    fx = f(x);
    if x~=0, ea = abs((x-x0)/x); end
    iter = iter + 1;
    results(iter+1,:) = [iter,x,fx,ea];
    
    % convert symbolic to double
    x = vpa(x,16); fx = vpa(fx,16);
    x0 = x; f0 = fx;
end

if iter < iter_max, results = results(1:iter+1,:); end

%% Printout

if ~silent
    disp(' ')
    disp('----------------------- Newton-Raphson ------------------------')
    disp(['Equation             : ',char(f)])
    disp(['Starting Point       : ',num2str(x00)])
    disp(['Optimal Result       : ',num2str(double(x))])
    disp(['Final/Expected Error : ',num2str(double(ea)*100,3),'%','/',num2str(es*100,3),'%'])
    disp(['Final/Max Iteration  : ',num2str(iter),'/',num2str(iter_max)])
    disp('Iterative Process    :')
    disp(array2table(results,'VariableNames',{'iter','x','fx','ea'}))
    disp('---------------------------------------------------------------')
end