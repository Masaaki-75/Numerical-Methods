function val = MidpointInt(x,f,npoint)

lb = min(x); ub = max(x);
[midpoints,h] = Midpoints(lb,ub,npoint);
val = sum(f(midpoints))*h;
