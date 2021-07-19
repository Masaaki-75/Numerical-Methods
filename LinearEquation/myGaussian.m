function [x,Ab] = myGaussian(varargin)
%==========================================================================
%MYGAUSSIAN Solving linear equation system using Gaussian method.
%
% Author: Ma Chenglong (github: Masaaki-75)
%         Copyright 2021 Sun Yat-sen University.
%
% Syntax: 
%
%   [X,AB] = MYGAUSSIAN(AUG);
%
%       returns the solution X as a column vector and the final augmented
%       matrix AB (could be simplest row echelon or row echelon) according
%       to the input augumented matrix AUG.
%
%   [X,AB] = MYGAUSSIAN(A,B);
%
%       returns the same results as above but takes coefficient matrix A
%       and constant column vector B as input.
%
%       If input argument B is a numeric scalar, then the function treats
%       it as a column vector whose elements take the same value of B.
%
%
%==========================================================================

%% Input Arguments Processing
Ab = varargin{1};
rows = size(Ab,1);

if nargin == 2
    A = Ab; b = varargin{2};
    if isscalar(b), b = b*ones(rows,1); end
    if size(b,2)>size(b,1), b = b.'; end
    Ab = [A,b];
end

%% Computation

% Judge the existence of solution.
A = Ab(:,1:end-1);
rankA = rank(A);
rankAb = rank(Ab);
disp(' ');

if rankA ~= rankAb  % case 1: NO solution
    x = A\b;
    disp('The input linear equation system has no solution.')
    return
else
    n_var = size(A,2);
    if rankA < n_var  % case 2: infinite solution
        disp('The input linear equation system has infinite number of solution.')
        disp('Particular solution and general solution are given respectively')
        disp('as two column vectors of the output array.')
        Ab = RowSimplest(Ab);  % simplest row echelon form
        x_p = Ab(:,end);  % particular solution
        x_g = null(Ab(:,1:end-1),'r');  % general solution
        x = [x_p, x_g];
        
    elseif rank(A) == n_var  % case 3: one solution
        disp('The input linear equation system has a unique solution.')
        % Apply Gaussian elimination to obtain row echelon form
        Ab = Eliminate(Ab);
        
        if n_var < rows  % case 3.1: redundant equation 
            Ab = ReduceRR(Ab);  % row simplest form
            x = Ab(1:end-1,end);
        elseif n_var == rows  % case 3.2: just right
            % backward substitution
            x = Substitute(Ab(:,1:end-1),Ab(:,end),'upper');
        end
    end
end