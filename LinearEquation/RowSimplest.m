function [Ar,ranks] = RowSimplest(A)
%==========================================================================
% Obtain the simplest (reduced) row echelon form of the input matrix.
%
% Author: Ma Chenglong (github: Masaaki-75)
%         Copyright 2021 Sun Yat-sen University.
%
% Syntax: 
%
%   [AR,RANKS] = ROWSIMPLEST(A);
%
%       returns the simplest row echelon form AR and the rank RANKS of the
%       input matrix A.
%
%==========================================================================
Ar = ReduceRR(Eliminate(A));

ranks = sum(any(Ar,2));
end
