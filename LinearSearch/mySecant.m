function [x,results] = mySecant(f,init,es,iter_max,modified,silent)
%==========================================================================
%MYSECANT estimate the root of univariate equation f(x)=0 using
% secant method.
%
% Author: Ma Chenglong (github: Masaaki-75)
%         Copyright 2021 Sun Yat-sen University.
%
% Syntax:
%
%   [X,RESULTS] = MYSECANT(F,INIT[,ES,ITER_MAX,MODIFIED,SILENT]);
%
%       returns the estimated zero X and a matrix RESULTS storing the
%       iteration results using method specified by MODIFIED, given input 
%       function F, initial condition INIT, error tolerance ES and maximum
%       number of iteration ITER_MAX. The last argument SILENT determines 
%       whether to display the iteration results in command line window.
%
%       F is the given function, whose class can be one of the following
%       {'string','char','function_handle','sym','symfun'}. This argument
%       will be treated as a symbolic function in MYNEWTON() function.
%
%       INIT is a numeric vector for initial condition. If MODIFIED ==
%       true, then INIT should include an initial point X00 and a positive
%       scalar for perturbation D, i.e.
%
%           INIT = [X00, D];
%
%       Otherwise, INIT should include two initial points X00 and X01, i.e.
%
%           INIT = [X00,X01];
%
%       ES is the preset absolute relative error as error tolerance.
%       Optional, 1e-5 as default.
%
%       ITER_MAX is a positive integer representing the maximum number of
%       iteration. Optional, 200 as default.
%
%       MODIFIED is a logical or char value determining whether to use
%       modified version of secant method. 
%
%       If MODIFIED == true, then the iteration will be formulated as
%
%           X(i+1) = X(i) - D*X(i)*F(X(i))/(F(X(i)+D*X(i))-F(X(i)));
%
%       with D the perturbation in X(i). Otherwise, the classical secant
%       method will be applied (as default case):
%
%           X(i+1) = X(i) - (X(i)-X(i-1))*F(X(i))/(F(X(i))-F(X(i-1)));
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
if ~exist('es','var')||isempty(es), es = 1e-5; end
if ~exist('iter_max','var')||isempty(iter_max), iter_max = 200; end
if ~exist('modified','var')||isempty(modified), modified = ' '; end
if ~exist('silent','var')||isempty(silent), silent = false; end
silent = lower(convertStringsToChars(silent));
modified = lower(convertStringsToChars(modified));
if ~islogical(silent)
    switch silent
        case {'on','disp','display','show','false'}
            silent = false;
        otherwise
            silent = true;
    end
end

%% Initial Setup
iter = 1;
% The "0-th" iteration for secant method is somehow hard to define,
% since we simply use two initial points to compute the first new point,
% which is recognized as the 1st iteration. 
% So here ITER = 1 (instead of ITER = 0). Accordingly, the RESULTS will
% take up ITER_MAX rows (instead of ITER_MAX + 1).
results = zeros(iter_max,6);
ea = 1;

%% Iteration

switch modified
    
    % Modified version of secant method
    case {'modified','delta','advanced','perturbation',1,true}
        
        modified = true;
        x00 = init(1); delta = init(2);
        
        if x00 == 0
            x00 = delta;  % x00 should be set non-zero
            warning(['Warning: Initial point should not be 0 in MODIFIED ',...
                'secant method. New initial point reset as ',num2str(x00),'.'])
        end
        
        % 1st iteration has begun (x is a symbolic variable):
        x = x00 - delta*x00*f(x00)/(f(x00+delta*x00)-f(x00));
        if isinf(x)||isnan(x), error('Error: Division by 0 when updating X.'); end
        % isa() function also works for symbolic variable
        results(1,:) = [iter,x00,delta,x,f(x),ea];  % store the 1st iteration
        
        while (ea > es)&&(iter < iter_max)
            iter = iter + 1;
            x00 = x;
            x = x00 - delta*x00*f(x00)/(f(x00+delta*x00)-f(x00));
            ea = abs((x-x00)/x);
            if isinf(ea)||isnan(ea), error('Error: Division by 0 when updating EA.'); end
            
            results(iter,:) = [iter,x00,delta,x,f(x),ea];
            
            x = vpa(x,16);  % retain the precision
            ea = double(ea);  % convert to double for logical calculation
        end
        
        
    % Vanilla version of secant method    
    otherwise  
        
        modified = false;
        x00 = init(1); x01 = init(2);
        
        % 1st iteration has begun (x is a symbolic variable):
        x = x01 - f(x01)*(x00-x01)/(f(x00) - f(x01));
        if isinf(x)||isnan(x), error('Error: Division by 0 when updating X.'); end
        
        results(1,:) = [iter,x00,x01,x,f(x),ea];   % store the 1st iteration
        
        while (ea > es) && (iter < iter_max)
            iter = iter + 1;
            x00 = x01; x01 = x;
            x = x01 - f(x01)*(x00-x01)/(f(x00) - f(x01));
            ea = abs((x-x01)/x);
            if isinf(ea)||isnan(ea), error('Error: Division by 0 when updating EA.'); end
            
            results(iter,:) = [iter,x00,x01,double(x),f(double(x)),ea];
            
            x = vpa(x,16); x01 = vpa(x01,16);   % retain the precision
            ea = double(ea);
        end
end

if iter<iter_max, results = results(1:iter,:); end  % clear redundant rows

%% Printout

if ~silent
    disp(' ')
    disp('--------------------------- Secant ----------------------------')
    disp(['Function             : ',char(f)])
    if modified
        disp(['Starting Point       : ',num2str(init(1))])
        disp(['Perturbation         : ',num2str(init(2))])
        temp1 = 'x0'; temp2 = 'delta';
    else
        disp(['Starting Point 1     : ',num2str(init(1))])
        disp(['Starting Point 2     : ',num2str(init(2))])
        temp1 = 'x00'; temp2 = 'x01';
    end
    disp(['Optimal Result       : ',num2str(double(x))])
    disp(['Final/Expected Error : ',num2str(ea*100,3),'%','/',num2str(es*100,3),'%'])
    disp(['Final/Max Iteration  : ',num2str(iter),'/',num2str(iter_max)])
    disp('Iterative Process    :')
    disp(array2table(results,'VariableNames',...
        {'iter',temp1,temp2,'x','fx','ea'}));
    disp('---------------------------------------------------------------')
end
