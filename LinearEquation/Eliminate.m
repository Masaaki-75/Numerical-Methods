function [Ae,row_ch_count] = Eliminate(A)
%==========================================================================
%ELIMINATE Row echelon form of the input matrix.
%
% Author: Ma Chenglong (github: Masaaki-75)
%         Copyright 2021 Sun Yat-sen University.
%
% Syntax: 
%
%   [AE,RCH] = ELIMINATE(A);
%
%       returns the row echelon form AE and the number of times of row
%       interchange caused by partial pivoting with input matrix A.
%
%==========================================================================

rows = size(A,1);
row_ch_count = 0;

for ii = 1 : rows - 1
    [A,temp] = PartialPivot(A,ii);  % partial pivoting
    row_ch_count = row_ch_count + temp;
    for jj = ii + 1 : rows
%         if abs(Aug(ii,ii)) < tol
%             Aug([ii,ii+1],:) = Aug([ii+1,ii],:);  % pivot interchanging
%             row_change_count = row_change_count + 1;
%         end
        facto = A(jj,ii)/A(ii,ii);
        % if "0/0" exists, the computation turns to next column.
        if isnan(facto), facto = A(jj,ii+1)/A(ii,ii+1); end
        A(jj,ii:end) = A(jj,ii:end) - A(ii,ii:end) * facto;  % elimination
    end
end

Ae = A;
Ae(abs(Ae)<=1e-10) = 0;  % remove the perturbation.
end