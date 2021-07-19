function metric = FittingError(y_hat,y,dof)
%==========================================================================
%FITTINGERROR evaluation of fitting results.
%
% Author: Ma Chenglong (github: Masaaki-75)
%         Copyright 2021 Sun Yat-sen University.
%
% Syntax: 
%
%   METRIC = FITTINGERROR(Y_HAT,Y,DOF);
%
%       returns evaluation results of the predicted values Y_HAT and the
%       obeserved values Y, with given freedom of degree DOF.
%
%==========================================================================

if size(y_hat,1) ~= size(y,1), y_hat = y_hat.'; end

n = length(y); assert(n==length(y_hat));

RSS = sum(norm((y_hat - y),2));
ESS = sum(norm(y_hat - mean(y),2));
r2 = 1-RSS/ESS;
stderr = sqrt(RSS/dof);
metric = [RSS;ESS;stderr;r2];