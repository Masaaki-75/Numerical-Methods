function [x,ea,iter] = myFPiter(f,x0,es,iter_max)

if ~exist('es','var'), es = 1e-5; end
if ~exist('iter_max','var'), iter_max = 200; end

iter = 0;
x = f(x0) + x0;
ea = 1;

while (ea>es)&&(iter<=iter_max)
    x0 = x;
    iter = iter + 1;
    x = f(x0) + x0;
    ea = abs((x-x0)/x);
    if isinf(ea)||isnan(ea), error('Error: Division by 0 from EA.'); end
end
