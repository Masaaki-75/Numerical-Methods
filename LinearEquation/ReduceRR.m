function Ar = ReduceRR(Ae)
%==========================================================================
%REDUCERR Simplest (reduced) row echelon form of the input row echelon
% matrix. 
%
% Note that this function serves as a intermediate part of ROWSIMPLEST()
% function, so the input matrix should ALREADY be a row echelon.
%
% Author: Ma Chenglong (github: Masaaki-75)
%         Copyright 2021 Sun Yat-sen University.
%
% Syntax: 
%
%   AR = REDUCERR(AE);
%
%       returns the simplest row echelon form AR of the input row echelon
%       matrix AE.
%
%==========================================================================

rows = size(Ae,1);
Ar = Ae;
Ar(abs(Ar)<=1e-10) = 0;  % remove perturbation.

% find pivot position.
loc = zeros(1,rows);
for ii = 1:rows
    temp = find(Ar(ii,:)~=0,1);
    if ~isempty(temp), loc(ii) = temp; end
end

% LOC stores the index of the pivot in each row.
% "LOC(ii) = 0" indicates no pivot in ii-th row.

for ii = rows:-1:2
    if loc(ii) ~= 0
        Ar(ii,loc(ii):end) = Ar(ii,loc(ii):end)/Ar(ii,loc(ii));  % pivot=1
        for jj = ii-1:-1:1
            facto = Ar(jj,loc(ii))*Ar(ii,loc(jj):end);
            Ar(jj,loc(jj):end) = Ar(jj,loc(jj):end) - facto;
        end
    end
end

if loc(1)~=0, Ar(1,:) = Ar(1,:)/Ar(1,loc(1)); end

Ar(abs(Ar)<=1e-10) = 0;
end