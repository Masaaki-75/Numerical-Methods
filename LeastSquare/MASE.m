function err = MASE(x,y,dim,type)
%==========================================================================
%MASE mean absolute/square error of the input data.
%
% Author: Ma Chenglong (github: Masaaki-75)
%         Copyright 2021 Sun Yat-sen University.
%
% Syntax: 
%
%   ERR = MASE(X,Y[,TYPE]);
%
%       returns the mean absolute error MAE or mean square error MSE of the
%       input data X and Y.
%
%       DIM is a positive integer specifying the dimension along which 
%       the mean value will be taken.
%
%       TYPE is a char specifying error metric. Optional. 'mse' as default.
%
%       - 'mse': mean square error.
%
%       - 'mae': mean absolute error.
%
%==========================================================================

% Default value setting
if size(x,1) ~= size(y,1), x = x.'; end
n = length(y); assert(n==length(x));

if ~exist('dim','var'), dim = 1; end 
if ~exist('type','var'), type = 'mse'; end

switch type
    case {'mse','MSE','square'}
        err = mean((x - y).^2,dim);
    case {'mae','MAE','abs'}
        err = mean(abs(x - y),dim);
    otherwise
        error('InputError: Argument "type" should either be "mse" or "mae"');
end