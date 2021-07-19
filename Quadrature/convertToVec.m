function out = convertToVec(in,dim)

siz = size(in);
n = numel(in);

if ~ismatrix(in)
    if isnumeric(in)  % input multi-channel tensor
        in = reshape(in,[1,n]);
    else
        out = in; return;
    end
end

if ~exist('dim','var')||isempty(dim), dim = 'row'; end
dim = convertStringsToChars(dim);

if isnumeric(dim)
    if dim == 1, dim = 'col';
    else, dim = 'row';
    end
elseif ischar(dim)
    dim = lower(dim);
end

dim = dim(1);

if strcmp(dim,'r')
    if min(siz) == 1
        out = reshape(in,[1,n]);
    else
        out = reshape(in.',[1,n]);
    end
else
    out = reshape(in,[n,1]);
end
