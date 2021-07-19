function [xopt,results] = myGSM(f,init,es,iter_max,silent)
%==========================================================================
%MYGSM estimate the critical point of univariate function f(x)
% using golden section method.
%
% Author: Ma Chenglong (github: Masaaki-75)
%         Copyright 2021 Sun Yat-sen University.
%
% Syntax:
%
%   [X,RESULTS] = MYGSM(F,INIT[,ES,ITER_MAX,SILENT]);
%
%       returns the estimated crtical point X and a matrix RESULTS storing
%       the relevant parameters computed in each iteration, given input 
%       function F, initial domain INIT, error tolerance ES and maximum
%       number of iteration ITER_MAX. The last argument SILENT determines 
%       whether to display the iteration results in command line window.
%
%       F is the given function, whose class can be one of the following
%       {'string','char','function_handle','sym','symfun'}. This argument
%       will be treated as a function handle in MYGSM() function.
%
%       INIT is a 2D vector specifying the initial domain (lower and upper
%       bounds included).
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

%% Initial Setup
R = 0.618033988749895;
init = sort(init);
xlow0 = init(1); xup0 = init(2);

%% Iteration
% first iteration.
xlow = xlow0; xup = xup0;
d = xup - xlow;
x1 = xlow + R*d; f1 = f(x1);
x2 = xup - R*d; f2 = f(x2);

if f1 > f2, xopt = x1; fopt = f1;
else, xopt = x2; fopt = f2; end

if xopt ~= 0, ea = (1-R)*abs((xup-xlow)/xopt); else, ea = 1; end

results = zeros(iter_max,11);
results(1,:) = [xlow,f(xlow),xup,f(xup),x2,f2,x1,f1,xopt,fopt,ea];

iter = 1;

while (ea > es) && (iter < iter_max)
    iter = iter + 1;
    if f1 > f2
        xlow = x2;  % update x_low
        d = xup - xlow;  % update stride
        xopt = x1; fopt = f1;  % update x_opt
        
        x2 = x1; f2 = f1;  % update x2_new = x1_old
        x1 = xlow + R*d;  % update x1_new
        f1 = f(x1);  % update f(x1)_new
    else
        xup = x1;  % update x_up
        d = xup - xlow;  % update stride
        xopt = x2; fopt = f2;  % update x_opt
        
        x1 = x2; f1 = f2;
        x2 = xup - R*d;
        f2 = f(x2);
    end
    
    if xopt ~= 0, ea = (1-R)*abs((xup-xlow)/xopt); end
    
    results(iter,:) = [xlow,f(xlow),xup,f(xup),x2,f2,x1,f1,xopt,fopt,ea];
end

if iter < iter_max, results = results(1:iter,:); end


%% Printout

if ~silent
    disp(' ')
    disp('----------------------- Golden Section ------------------------')
    disp(['Function             : ',char(f)])
    disp(['Initial Scope        : ','[',num2str(xlow0),',',num2str(xup0),']'])
    disp(['Final Scope          : ','[',num2str(xlow),',',num2str(xup),']'])
    disp(['Optimal Result       : ',num2str(xopt)])
    disp(['Final/Expected Error : ',num2str(ea*100,3),'%','/',num2str(es*100,3),'%'])
    disp(['Final/Max Iteration  : ',num2str(iter),'/',num2str(iter_max)])
    disp('Iterative Process    :')
    disp(array2table(results,'VariableNames',...
        {'xl','fl','xu','fu','x2','f2','x1','f1','xop','fop','ea'}));
    disp('----------------------------------------------------')
end
