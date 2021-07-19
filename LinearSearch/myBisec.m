function [x,results] = myBisec(f,init,es,iter_max,silent)
%==========================================================================
%MYBISEC estimate the root of univariate equation f(x)=0 using
% bisection methods.
%
% Author: Ma Chenglong (github: Masaaki-75)
%         Copyright 2021 Sun Yat-sen University.
%
% Syntax:
%   [x[,err,iter]] = myBisec(f,xl,xu[,es,iter_max]);
%
% Input:
% - f       : Given function (@(x)(func) as input).
%
% - xl      : Lower bound estimation of the root.
%
% - xu      : Upper bound estimation of the root.
%
% - es      : Preset absolute relative error as error tolerance.
%             OPTIONAL, 5% as default.
%
% - iter_max: Maximum numbers of iteration.
%             OPTIONAL, log2[(xu-xl)/es] as default.
%
% Output:
% - x       : Estimated root.
%
% - err     : Final absolute relative error.
%
% - iter    : Final iteration numbers.
%==========================================================================


%% Input Arguments Processing
f = convertStringsToChars(f);
if isa(f,'char'), f = str2sym(f); end
if isa(f,'sym')||isa(f,'symfun'), f = matlabFunction(f); end
if ~exist('es','var')||isempty(es), es = 1e-5; end
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
init = sort(init);
if ~exist('iter_max','var') || isempty(iter_max)
    iter_max = ceil(log2((init(2)-init(1))/es))+1;
end

%% Initial Setup
xl = init(1); xu = init(2);
iter = 0;  % initial val of current # of iteration
x_old = xl;  % initial val of x_old
ea = 1;  % initial val of ARE
fl = f(xl);
results = zeros(iter_max+1,8);
results(1,:) = [iter,xl,fl,xu,f(xu),x_old,fl,ea];

%% Iteration
while (ea > es) && (iter < iter_max)
    x = (xl + xu)/2;
    f_new = f(x);
    iter = iter + 1;
    if x ~= 0
        ea = abs((x - x_old)/x);
    end
    
    % Update margin [xl, xu]
    crit = fl*f_new;
    if crit < 0  % root exists in [xl, x_new]
        x_old = x;
        xu = x;  % update upper bound
    elseif crit > 0  % root exists in [x_new, xu]
        x_old = x;
        xl = x;  % update lower bound
        fl = f_new;
    else  % x_new = root
        ea = 0;
    end
    
    results(iter+1,:) = [iter,xl,fl,xu,f(xu),x,f_new,ea];
end

if iter < iter_max, results = results(1:iter+1,:); end

%% Printout
if ~silent
    disp(' ')
    disp('-------------------------- Bisection --------------------------')
    disp(['Function             : ',char(f)])
    disp(['Initial Scope        : ','[',num2str(init(1)),',',num2str(init(2)),']'])
    disp(['Final Scope          : ','[',num2str(xl),',',num2str(xu),']'])
    disp(['Optimal Result       : ',num2str(x)])
    disp(['Final/Expected Error : ',num2str(ea*100,3),'%','/',num2str(es*100,3),'%'])
    disp(['Final/Max Iteration  : ',num2str(iter),'/',num2str(iter_max)])
    disp('Iterative Process    :')
    disp(array2table(results,'VariableNames',{'iter','xl','fl','xu','fu','x','fx','ea'}))
    disp('---------------------------------------------------------------')
end