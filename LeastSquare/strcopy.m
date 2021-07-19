function str_new = strcopy(varargin)
%==========================================================================
%STRCOPY replicate the input string for specific times.
%
% Author: Ma Chenglong (github: Masaaki-75)
%         Copyright 2021 Sun Yat-sen University.
%
% Syntax: 
%
%  str_new = strcopy(str,n);
%
%       returns a n-time replication of the input string. For example,
%           strcopy('ab-',3)
%       returns 'ab-ab-ab-'.
%
%==========================================================================
temp1 = varargin{1};
temp2 = varargin{2};

if ischar(temp1)
    validateattributes(temp2,{'numeric'},{'scalar'});
    str = temp1; n = temp2;
elseif isnumeric(temp1)
    validateattributes(temp2,{'char'},{'scalartext'});
    str = temp2; n = temp1;
end

str_new = str;

for ii = 1:n-1
    str_new = strcat(str_new,str);
end