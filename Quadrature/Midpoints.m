function [midpoints,step] = Midpoints(lb,ub,npoint)

step = (ub-lb)/npoint;
m1 = lb + step/2;
m2 = ub - step/2;
midpoints = linspace(m1,m2,npoint);
