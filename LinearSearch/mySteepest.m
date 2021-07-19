function [crit,results] = mySteepest(f,init,varargin)
%==========================================================================
%MYSTEEPST estimate the critical point where the input function takes
% its extremum.
%
% Author: Ma Chenglong (github: Masaaki-75)
%         Copyright 2021 Sun Yat-sen University.
%
% Syntax:
%
%   [CRIT,RESULTS] = MYSTEEPEST(F,INIT,...
%                               'MinGradNorm',GNORM_MIN,...
%                               'MaxIter',ITER_MAX,...
%                               'Method',METHOD,...
%                               'Param',PARAM,...
%                               'Display',DISPTYPE);
%
%       returns the critical point CRIT and iterative results as a matrix
%       RESULTS.
%
%       F is the given function, whose class can be one of the following
%       {'string','char','function_handle','sym','symfun'}.
%
%       INIT is a 2D vector specifying the initial point.
%
%       Other arguments are name-value pair arguments:
%
%       GNORM_MIN is a numeric positive scalar representing the minimum
%       value of the norm of the gradient. Optional, 1e-5 as default.
%
%       ITER_MAX is a positive integer representing the maximum number of
%       iteration. Optional, 400 as default.
%
%       METHOD is string or char specifying the search strategy of optimal
%       stride. Possible values are listed below. 
%       Optional, 'none' as default.
%
%       - 'none': fixed stride. Here, the argument PARAM should be a 1x1
%         cell containing the constant stride.
%
%       - 'rate': the value of stride decreases as the gradient norm
%         decreases.
%
%       - 'newton': use Newton-Raphson algorithm to search the optimal
%         stride. Here, PARAM should contain the following
%         arguments:
%         {initial point, error tolerance, max # of iteration}.
%
%       - 'gold': use golden section algorithm to search the optimal
%         stride. Here, PARAM should contain the following
%         arguments:
%         {initial domain, error tolerance, max # of iteration}.
%
%       PARAM is a cell containing necessary arguments of the selected
%       seach strategy and is dependent on METHOD.
%
%       DISPTYPE is a char specifying how the iteration results are
%       presented. Possible values are listed below. 
%       Optional, 'figure' as default.
%
%       - 'figure': Show the contour and surface map of the function and
%         along with the iterative process.
%
%       - 'table': Show the iteration results in the form of a table in the
%         command line window.
%
%       - 'both': Show figure and table simultaneously.
%        
%       - 'none': Do not show the results.
%
%==========================================================================

%% Input Arguments Processing
defaultMethod = 'rate';
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
validFunction = @(x) isa(x,'char')||isa(x,'function_handle')...
    ||isa(x,'sym')||isa(x,'symfun');

p = inputParser;
addRequired(p,'Function', validFunction);  % required argument
addRequired(p,'InitPoint', @(x) isnumeric(x));
addOptional(p,'MinGradNorm',1e-5,validScalarPosNum);  % optional argument
addOptional(p,'MaxIter',200,validScalarPosNum);
addParameter(p,'Method',defaultMethod,...
    @(x) ischar(convertStringsToChars(x)));
addParameter(p,'Param',[],@(x) iscell(x));
addParameter(p,'Display','figure',@(x) ischar(x)||isstring(x));
parse(p,f,init,varargin{:});

f = p.Results.Function;
init = p.Results.InitPoint;
gnorm_min = p.Results.MinGradNorm;
iter_max = p.Results.MaxIter;
method = p.Results.Method;
param = p.Results.Param;
disptype = p.Results.Display;

% Default values setting
if isa(f,'function_handle'), f = sym(f); end
if isa(f,'char'), f = str2sym(f); end
if isa(f,'sym'), f = symfun(f,symvar(f)); end
if ~exist('gnorm_min','var')||isempty(gnorm_min), gnorm_min = 1e-5; end
if ~exist('iter_max','var')||isempty(iter_max), iter_max = 400; end
if ~exist('method','var')||isempty(method)||~exist('param','var')||isempty(param)
    method = 'rate'; param = {1};
end
if ~exist('disptype','var')||isempty(disptype), disptype = 'figure'; end
method = lower(convertStringsToChars(method));
disptype = lower(convertStringsToChars(disptype));

%% Initial Setup
iter = 0;  % count of iteraion
results = zeros(iter_max,7);  % pre-allocate memory to store results.
x1 = init(1); y1 = init(2);  % x and y of the initial point
vars = symvar(f);  % variables of the given function
x = vars(1); y = vars(2);
h = sym('h','real');  % stride

fx = diff(f,x); fy = diff(f,y);  % partial derivative
% compute gradient value.
px = double(fx(x1,y1));
py = double(fy(x1,y1));
gnorm = sqrt(px.*px+py.*py);  % Euclidean norm of gradient

if gnorm <= gnorm_min, return; end

%% Iteration
while gnorm > gnorm_min && iter < iter_max
    iter = iter + 1;
    gx = px./gnorm; gy = py./gnorm;
    x1_s = x1 - gx * h;
    y1_s = y1 - gy * h;
    
    switch method
        case {'newton'}  % Newton-Raphson method
            fh = f(x1_s,y1_s);
            % Here x1_s, y1_s & fh are all symbolic variables.
            % But when fh is entered into myNewton() function, it would be
            % converted into a symbolic function.
            
            h_opt = myNewton(diff(fh),param{1},param{2},param{3},param{4});
            % h_opt is a double variable.
        case {'gold'}  % golden section method
            fh = f(x1_s,y1_s);
            h_opt = myGSM(-fh,param{1},param{2},param{3},param{4});
        case {'rate'}  % moving stride method
            
            if gnorm <= 1e-1 && gnorm > 1e-2, lr = 1e-1;
            elseif gnorm <= 1e-2 && gnorm > 1e-3, lr = 1e-2;
            elseif gnorm <= 1e-3 && gnorm > 1e-4, lr = 1e-3;
            elseif gnorm <= 1e-4 && gnorm > 1e-5, lr = 1e-4;
            elseif gnorm <= 1e-5, lr = 1e-5;
            else, lr = 1;
            end
            h_opt = param{1}*lr;
            
        otherwise
            h_opt = param{1};
    end
    
    % update critical point.
    x1 = double(subs(x1_s,h,h_opt));
    y1 = double(subs(y1_s,h,h_opt));
    
    % update gradient.
    px = double(fx(x1,y1));
    py = double(fy(x1,y1));
    gnorm = sqrt(px.*px+py.*py);
    
    % store the data.
    results(iter,:) = [iter,h_opt,x1,y1,px,py,gnorm];
    
    
end

if iter < iter_max, results = results(1:iter,:); end
crit = [x1,y1];

%% Visualization
switch disptype
    case {'figure','fig','plot'}
        showFigure(f,init,crit,results);
    case {'text','txt','table','tab','print'}
        showTable(f,init,crit,gnorm,gnorm_min,iter,iter_max,results);
    case {'all','show','both'}
        showTable(f,init,crit,gnorm,gnorm_min,iter,iter_max,results);
        showFigure(f,init,crit,results);
    otherwise
        disp(' ')
end

end

function showFigure(f,init,crit,results)
x1 = crit(1); y1 = crit(2);
x0 = init(1); y0 = init(2);
eqn = char(f); eqn = replace(eqn,'exp','{\rm exp}');
eqn = replace(eqn,'log','{\rm log}');
eqn = replace(eqn,'*','');

figure,
xmin = min(x1,x0); xmax = max(x1,x0);
ymin = min(y1,y0); ymax = max(y1,y0);
eqntemp = replace(eqn,'x','x_1');
eqntemp = replace(eqntemp,'y','x_2');

subplot(1,2,1)
fcontour(f,[xmin-2 xmax+2 ymin-2 ymax+2],'linewidth',2),
hold on, grid on
xlabel('$x_1$','interpreter','latex','fontsize',14)
ylabel('$x_2$','interpreter','latex','fontsize',14)
title(['Contour of $f({\rm \bf x})=',eqntemp,'$'],'interpreter','latex','fontsize',14)
%quiver(init(1),init(2),1.2*init(1),1.2*init(2));
%plot(init(1),init(2),'*','color',[0.5,0,0],'markersize',2)
points = [x0,y0];
points = [points;[results(:,3),results(:,4)]];
np = size(points,1);
colors = [ones(np,1),linspace(0,1,np)',zeros(np,1)];

subplot(1,2,2)
fsurf(f,[xmin-2 xmax+2 ymin-2 ymax+2],'EdgeColor','none'), hold on
xlabel('$x_1$','interpreter','latex','fontsize',14)
ylabel('$x_2$','interpreter','latex','fontsize',14)
zlabel('$f({\bf \rm x})$','interpreter','latex','fontsize',14)
title(['Surface of $f({\rm \bf x})=',eqntemp,'$'],'interpreter','latex','fontsize',14)
for ii = 1:np-1
    subplot(1,2,1)
    plot([points(ii,1),points(ii+1,1)],[points(ii,2),points(ii+1,2)],...
        '*--','color',colors(ii,:),'LineWidth',1.2,'MarkerSize',5)
    %drawarrow(points(ii,:),points(ii+1,:),'arrow',[red(ii),0,1])
    subplot(1,2,2)
    val1 = double(f(points(ii,1),points(ii,2)));
    val2 = double(f(points(ii+1,1),points(ii+1,2)));
    plot3([points(ii,1),points(ii+1,1)],[points(ii,2),points(ii+1,2)],[val1,val2],...
        'o--','color',colors(ii,:),'LineWidth',1.2,'MarkerSize',5,...
        'MarkerFaceColor',colors(ii,:),'MarkerEdgeColor',colors(ii,:))
end
end

function showTable(f,init,crit,gnorm,gnorm_min,iter,iter_max,results)
disp(' ')
disp('------------------ Steepest Gradient Descent ------------------')
eqn = char(f); eqn = replace(eqn,'*','\times ');
disp(['Function                     : ',eqn])
disp(['Starting Point [x0,y0]       : ','[',num2str(init(1)),',',num2str(init(2)),']'])
disp(['Optimal Result [x1,y1]       : ','[',num2str(crit(1)),',',num2str(crit(2)),']'])
disp(['Final/Expected Gradient Norm : ',num2str(gnorm,3),'%','/',num2str(gnorm_min,3),'%'])
disp(['Final/Max Iteration          : ',num2str(iter),'/',num2str(iter_max)])
disp('Iterative Process            :')
disp(array2table(results,'VariableNames',...
    {'iter','h','x','y','dfdx','dfdy','grad_norm'}));
disp('---------------------------------------------------------------')
end