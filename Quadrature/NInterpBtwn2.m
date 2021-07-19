function varargout = NInterpBtwn2(varargin)

if nargin == 2
    if is2dvector(varargin{1}) && isscalar(varargin{2})
        y = varargin{1}; N = varargin{2};
    else
        y = varargin{2}; N = varargin{1};
    end
elseif nargin == 3
    if is2dvector(varargin{1}) && is2dvector(varargin{2})
        % input x,y and N
        x = varargin{1}; y = varargin{2}; N = varargin{3};
        if ~(isscalar(N) && isnumeric(N))
            error('Error: Number of interpolated points has not been specified.')
        end
    elseif is2dvector(varargin{1}) && ~is2dvector(varargin{2})
        % input y, N and method
        y = varargin{1};
        if isscalar(varargin{2}) && isnumeric(varargin{2})
            N = varargin{2}; method = varargin{3};
        else
            N = varargin{2}; method = varargin{3};
        end
    end
elseif nargin == 4
    % input x, y, N and method
    x = varargin{1}; y = varargin{2};
    if isscalar(varargin{3}) && isnumeric(varargin{3})
        N = varargin{3}; method = varargin{4};
    else
        N = varargin{4}; method = varargin{3};
    end
end

if ~exist('x','var')||isempty(x), x = (1:length(y))'; end
if ~exist('method','var')||isempty(method), method = 'linear'; end

x = convertToColVec(x);
y = convertToColVec(y);

n_in = length(x);
n_out = n_in + (n_in - 1) * N;

x1 = zeros(n_out,1);
y1 = zeros(n_out,1);

for ii = 1:n_in-1
    xi = x(ii:ii+1);
    yi = y(ii:ii+1);
    xq = linspace(xi(1),xi(2),N+2);
    xq = xq(1:end-1);
    yq = interp1(xi,yi,xq,method);
    
    idx1 = (ii-1)*N + ii;
    idx2 = ii*(N+1);
    x1(idx1:idx2) = xq;
    y1(idx1:idx2) = yq;
end

x1(end) = x(end);
y1(end) = y(end);


if nargout == 2
    varargout{1} = x1;
    varargout{2} = y1;
else
    varargout{1} = y1;
end

end

function vc = convertToColVec(v)

siz = size(v);

if siz(2) > siz(1) && siz(1) == 1
    vc = v.';
elseif min(siz(1),siz(2)) > 1
    n = numel(v);
    vc = reshape(v,n,1);
else
    vc = v;
end

end

function flag = is2dvector(a)

flag = false;
siz = size(a);

if isnumeric(a) && numel(siz)<=2 && min(siz)==1 && max(siz)>1
    flag = true;
end

end