clear; close all

f = '(x-3).^2 + (y-2).^2';
g = 'x*exp(-x^2-y^2)';
x0 = 1; y0 = 1;

gnorm_min = 1e-3;
[crit11,res11] = mySteepest(f,[x0,y0],...
    'MinGradNorm',gnorm_min,'MaxIter',100,'Method','Newton','Param',{1,1e-4,200,true},'Display','figure');
[crit12,res12] = mySteepest(f,[x0,y0],...
    'MinGradNorm',gnorm_min,'MaxIter',100,'Method','Gold','Param',{[0,2],1e-4,200,true},'Display','figure');
[crit13,res13] = mySteepest(f,[x0,y0],...
    'MinGradNorm',gnorm_min,'MaxIter',100,'Method','rate','Param',{0.15},'Display','figure');
[crit14,res14] = mySteepest(f,[x0,y0],...
    'MinGradNorm',gnorm_min,'MaxIter',450,'Method','none','Param',{0.0052},'Display','text');

x00 = 0; y00 = 0;
% if [x00,y00]==[1,1], the critical point turns out to fall far on the 
% plane where gradf(x,y) ~= 0 and ALSO f(x,y) ~= 0.
% This is not what we expected since we know that min{f(x,y)} < 0.
[crit21,res21] = mySteepest(g,[x00,y00],...
    'MinGradNorm',gnorm_min,'MaxIter',100,'Method','Newton','Param',{1,1e-4,200,true},'Display','figure');
[crit22,res22] = mySteepest(g,[.5,1],...
    'MinGradNorm',gnorm_min,'MaxIter',100,'Method','Gold','Param',{[0,2],1e-4,200,true},'Display','figure');
[crit23,res23] = mySteepest(g,[x00,y00],...
    'MinGradNorm',gnorm_min,'MaxIter',100,'Method','rate','Param',{0.1},'Display','figure');
[crit24,res24] = mySteepest(g,[x00,y00],...
    'MinGradNorm',gnorm_min,'MaxIter',100,'Method','none','Param',{0.01},'Display','text');