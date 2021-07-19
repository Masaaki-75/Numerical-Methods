function [Ach, row_ch_count] = PartialPivot(A,row)
%==========================================================================
%PARTIALPIVOT partial pivoting based on certain row of the matrix to avoid
% division by zero and reduce the round-off error.
%
% Author: Ma Chenglong (github: Masaaki-75)
%         Copyright 2021 Sun Yat-sen University.
%
% Syntax: 
%
%   [ACH, RCH] = PARTIALPIVOT(A,R)
%
%       returns the row echelon form AE and the number of times of row
%       interchange caused by partial pivoting with input matrix A.
%
%==========================================================================

row_ch_count = 0;
Ach = A;
rows = size(Ach,1);
row_current = row;
val_current = abs(Ach(row,row));

for ii = row + 1 : rows
    % search larger coeffiecient |Aii| in the following rows.
    temp = abs(Ach(ii,row));
    if temp > val_current
        val_current = temp;
        row_current = ii;
    end
end

if row_current ~= row
    % If larger |Aii| exists, interchange that row with current one.
    row_ch_count = row_ch_count + 1;
    Ach([row_current,row],row:end) = Ach([row,row_current],row:end);
end
end
    