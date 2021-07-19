clear; clc; close all
%% Data Preparation
x = (1900:10:2000).';
y = [75.995 91.972 105.711 123.203 131.669 150.697 179.323 203.212 226.505 249.633 281.422].';
x_new = [2010 2020].';
y_new = [308.745538 331.449281].';

%% Initial Setup
deg = 5;  % degree of polynomial
alpha = 1/100;  % scaling factor
scale = (alpha.^(-1:(deg-1))).';  % scaling vector

% unscaled augmented data
x_all = [x;x_new];
y_all = [y;y_new];

% scaled data ('p' stands for 'prime')
xp = x * alpha;  % Year 1990~2000 (every 10 yrs)
yp = y * alpha;
xp_all = x_all * alpha;  % Year 1990~2020 (every 10 yrs)
yp_all = y_all * alpha;

n = length(x);
n_new = length(x_new);
n_all = n + n_new;

% pre-allocating memory
ap = zeros(deg+1,deg);  % scaled coeffcient vector
yp_hat = zeros(n,deg);  % scaled predicted values (1990~2000)
metricp = ones(4,deg);  % fitting evaluation on scaled data (1990~2000)

a = zeros(deg+1,deg);  % unscaled coeffcient vector
y_hat_all = zeros(n_all,deg);  % unscaled predicted values (1990~2020)

y_hat_new = zeros(n_new,deg);  % unscaled predicted values (2010,2020)
err_new = ones(1,deg);  % errors evaluated on unscaled data (2010,2020)

%% Polynomial Regression
for m = 1:deg
    [ap(1:(m+1),m),yp_hat(:,m),metricp(:,m)] = myPolyReg(xp,yp,m,'fig');
    
    a(:,m) = ap(:,m).*scale;  % obtain unscaled coefficients
    y_hat_all(:,m) = polyval(flip(a(:,m)),x_all);
    % unscaled predicted values (1990~2020)
    
    y_hat_new(:,m) = polyval(flip(a(:,m)),x_new);
    % unscaled predicted values (2010,2020)
    
    err_new(m) = MASE(y_hat_new(:,m),y_new);
    
end

%% Visualization
figure  % show regression results in unscaled form
plot(x_all(1:end-2),y_all(1:end-2),'o','markerfacecolor','#0072BD',...
    'markeredgecolor','none','markersize',7)
hold on
plot(x_all(end-1:end),y_all(end-1:end),'o','markerfacecolor','r',...
    'markeredgecolor','none','markersize',7)
plot(x_all,y_hat_all,'linewidth',1)
xlabel('$t$/ year','interpreter','latex','fontsize',14)
ylabel('$y/\ \times 10^6$','interpreter','latex','fontsize',14)
title('Fitting curve of $y=\sum_{l=0}^{m}{a_l x^l}$',...
    'interpreter','latex','fontsize',14)
legend('raw data','augmented data','$m=1$','$m=2$','$m=3$','$m=4$','$m=5$',...
    'interpreter','latex','fontsize',14,'location','northwest')
