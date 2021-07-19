function [a,y_hat,metric,fhandle] = myPolyReg(x,y,deg,disptype)
%==========================================================================
%MYPOLYREG polynomial regression.
%
% Author: Ma Chenglong (github: Masaaki-75)
%         Copyright 2021 Sun Yat-sen University.
%
% Syntax: 
%
%   [A,Y_HAT,METRIC,FHANDLE] = MYPOLYREG(X,Y,DEG[,DISPTYPE]);
%
%       returns the coefficient vector A, predicted values Y_HAT,
%       evaluation results METRIC, polynomial function handle FHANDLE.    
%
%       X and Y are column vectors of data of independent and dependent
%       variable, respectively.
%
%       DEG is a positive integer representing the degree of polynomial.
%
%       DISPTYPE is a char specifying how the iteration results are
%       presented. Possible values are listed below. 
%       Optional, 'figure' as default.
%
%       - 'figure': Show the data points and the fitting curve.
%
%       - 'text': Print the fitting results in the command line window.
%
%       - 'both': Show figure and text simultaneously.
%        
%       - 'none': Do not show the results.
%
%==========================================================================

%% Input Arguments Processing
sizx = size(x); sizy = size(y);
if min(sizx) ~= 1 || min(sizy) ~= 1
    error('Error: Input data should be vectors.')
else
    if sizx(2) > sizx(1), x = x.'; end
    if sizy(2) > sizy(1), y = y.'; end
end
n = length(x); assert(n==length(y));
if ~exist('disptype','var')||isempty(disptype), disptype = 'none'; end
disptype = lower(convertStringsToChars(disptype));

%% Solving System of Linear Equations
X = Vandermonde(x,deg);
X1 = (X.')*X; Y1 = (X.')*y;

if abs(det(X1)) < 1e-10 || (rank(X1) ~= rank([X1,Y1]))
    warning("Warning: The matrix X'*X is singular (possibly non-invertible).")
    a = pinv(X1)*Y1;
else
    a = myGaussian(X1,Y1);
end

%% Error Evaluation
y_hat = X*a;
metric = FittingError(y_hat,y,n-deg-1);

%% Visualization
[fhandle,eqn] = PolyHandle(a);

switch disptype
    case {'text','txt','table','tab','print'}
        showText(deg,eqn,metric)
    case {'figue','fig','plot'}
        showFigure(x,y,y_hat,deg,metric)
    case {'all','show','both'}
        showText(deg,eqn,metric)
        showFigure(x,y,y_hat,deg,metric)
    otherwise
        disp(' ')
end

end


function showText(deg,eqn,metric)
disp(' ')
disp('-------------------- Polynomial Regression --------------------')
disp(['Degree of Polynomial       : ',num2str(deg)])
disp(['Fitting Curve Formula      : ',eqn])
disp(['Residual Sum Square (RSS)  : ',num2str(metric(1))])
disp(['Explained Sum Square (ESS) : ',num2str(metric(2))])
disp(['Standard Error of Estimate : ',num2str(metric(3))])
disp(['Determination Coefficient  : ',num2str(metric(4)*100),'%'])
disp('---------------------------------------------------------------')
end

function showFigure(x,y,y_hat,deg,metric)
figure
plot(x,y,'o','markerfacecolor','#0072BD'),hold on
plot(x,y_hat,'linewidth',1,'color','#D95319')
xlabel('$x$','interpreter','latex','fontsize',14)
ylabel('$y$','interpreter','latex','fontsize',14)
title(['Fitting curve of $y=\sum_{l=0}^{',num2str(deg),'}{a_l x^l}$'],...
    'interpreter','latex','fontsize',14)
legend('data points',['fitting curve',...
    newline,'($s_{y/x}$=',num2str(metric(3)),')',...
    newline,'($r^2=',num2str(metric(4)*100),'\%$)'],...
    'interpreter','latex','fontsize',14,'location','southeast')
end