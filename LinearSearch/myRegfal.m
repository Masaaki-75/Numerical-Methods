function [x,results] = myRegfal(f,init,es,iter_max,silent)
%==========================================================================
%MYREGFAL estimate the root of univariate equation f(x)=0 using
% Regula-Falsi method (or linear interpolation method).
%
% Author: Ma Chenglong (github: Masaaki-75)
%         Copyright 2021 Sun Yat-sen University.
%
% Syntax:
%   [x[,err,iter]] = myRegfal(f,xl,xu[,es,iter_max]);
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
%             OPTIONAL, 200 as default.
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
init = sort(init);

%% Initial Setup
xl = init(1); xu = init(2);
iter = 0; ea = 1;
x_old = xl;
fl = f(xl); fu = f(xu);
count_l = 0; count_u = 0;
results = zeros(iter_max+1,8);
results(1,:) = [iter,xl,fl,xu,fu,x_old,fl,ea];

%% Iteration
while (ea > es) && (iter < iter_max)
    
    % Halve the function value of stagnant boundary point
    % to boost the convergence. 
    if count_u >= 2, fu = fu/2; end
    if count_l >= 2, fl = fl/2; end
    
    x = xu - fu*(xu-xl)/(fu-fl);  % update X
    f_new = f(x);  % update function value
    iter = iter + 1;
    if x ~= 0
        ea = abs((x - x_old)/x);
    end
    
    crit = fl*f_new;
    
    % update margin [xl, xu]
    if crit < 0  % root exists in [xl, x_new]
        xu = x;  % update upper bound
        fu = f_new;
        x_old = x;
        count_l = count_l + 1;  % suspected stagnant point at xl
        count_u = 0;
        
    elseif crit > 0  % root exists in [x_new, xu]
        xl = x;  % update lower bound
        fl = f_new;
        x_old = x;
        count_u = count_u + 1;  % suspected stagnant point at xu
        count_l = 0;
        
    else  % x_new = root
        ea = 0;
    end
    
    results(iter+1,:) = [iter,xl,fl,xu,fu,x,f_new,ea];
end

if iter < iter_max, results = results(1:iter+1,:); end

%% Printout

if ~silent
    disp(' ')
    disp('------------------------ Regula Falsi -------------------------')
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