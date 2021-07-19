function [root,fac,remainder] = myBairstow(p,init,es,iter_max)
%==========================================================================
%MYBAIRSTOW estimates all the roots of univariate polynomial equation
% using Bairstow method. In other words, this function seeks to solve the
% following equation
%
%     F(X) = A(N)*X^N + A(N-1)*X^(N-1) + ... + A(1)*X + A(0) = 0
%
% by decomposing the N-degree polynomial F(X) into the product of quadratic
% factors Q(X) and a remainder polynomial R(X) with low degree as follows:
%
%     F(X) = Q(X)*R(X)
%     Q(X) = (X^2+R1*X+S1)*(X^2+R2*X+S2)*...*(X^2+Rm*X+Sm*X)
%     R(X) = some polynomial whose degree is <= 2.
%
% Author: Ma Chenglong (github: Masaaki-75)
%         Copyright 2021 Sun Yat-sen University.
%
% Syntax:
%
%   [R,FAC,RES] = MYBAIRSTOW(P,INIT[,ES,ITER_MAX]);
%
%       returns the estimated solutions (R), coefficient matrix (FAC) of
%       the decomposed quadratic factor FAC, and the coefficient vector
%       (RES) of the remainder polynomial, given polynomial equation
%       defined by coefficient vector P, initial guess (INIT) of
%       coefficients of quadratic factor, tolerance ES, maximum number of 
%       iteration ITER_MAX.
%
%       As described above, the coefficients of terms of Q(X) 
%           [[1,R1,S1];[1,R2,S2];...;[1,Rm,Sm]]
%       forms the output matrix FAC, and the output vector RES contains the
%       coefficients of R(X).
%
%       The input coefficient vector P is a row vector defining the
%       N-degree polynomial equation as described below:
%
%           P(1)*X^N + P(2)*X^(N-1) + ... + P(N)*X + P(N+1);
%
%       with the 1st entry containing the coefficient of the term with the
%       highest power.
%
%       INIT is a numeric vector for initial guess of coefficients of
%       quadratic factors, i.e. INIT = [R0,S0].
%
%       ES is the preset absolute relative error as error tolerance when
%       optimizing coefficients (R and S) of quadratic factor. Optional,
%       1e-5 as default.
%
%       ITER_MAX is a numeric vector representing the maximum number of
%       iteration. The 1st element is for the external iteration
%       (factorizing the polynomial) and the 2nd is for the internal
%       iteration (optimizing the quadratic coefficients). Optional,
%       [floor((N-1)/2),100] as default.
%
%==========================================================================

%% Input Arguments Processing
if ~exist('es','var')||isempty(es), es = 1e-5; end
if ~exist('init','var')||isempty(init), init = [-1,-1]; end
if min(size(p)) ~= 1
    error('Error: Input coefficient array should be a row vector.');
end
if size(p,1) ~= 1, p = p.'; end 

degree = size(p,2) - 1;  % degree of polynomial
if ~exist('iter_max','var')||isempty(iter_max)
    iter_max = [floor((degree-1)/2),100];
end
iter_fac_max = iter_max(1); iter_rs_max = iter_max(2);

%% Initial Setup
A = flip(p);
% Due to writing habits, input coefficients are usually entered as
%     p = [a(n),a(n-1),...,a(1),a(0)];
% However, in this function it would be better to arrange the coefficient
% vector in the following form:
%     A = [a(0),a(1),...,a(n-1),a(n)];
% hence the FLIP().

r = -init(1); s = -init(2);
% The calculation below treats the quadratic form as 
% (X^2 - R*X - S*X), so a minus sign is needed.

n = degree; iter_fac = 0;

% preallocating memory
root = zeros(degree,1);
fac = ones(iter_fac_max,3);


%% Iteration

while (n > 2) && (iter_fac < iter_fac_max)
    iter_fac = iter_fac + 1;
    B = zeros(1,n+1);  % same length as A
    C = zeros(1,n);
    
    ea_r = 1; ea_s = 1;
    % approxi. error of R and S
    iter_rs = 0;
    
    while (ea_r > es) || (ea_s > es) && (iter_rs < iter_rs_max)
        iter_rs = iter_rs + 1;
        
        B(n+1) = A(n+1);
        B(n) = A(n) + r * B(n+1);
        
        C(n) = B(n+1);
        C(n-1) = B(n) + r * C(n);
        
        for ii = n-1:-1:1
            B(ii) = A(ii) + r * B(ii+1) + s * B(ii+2);
            if ii <= n-2
                C(ii) = B(ii+1) + r * C(ii+1) + s * C(ii+2);
            end
        end
        
        deter = C(2)*C(2) - C(3)*C(1);
        if deter ~= 0
            dr = (B(1)*C(3) - B(2)*C(2))/deter;
            ds = (B(2)*C(1) - B(1)*C(2))/deter;
            r = r + dr; s = s + ds;
            if r ~= 0, ea_r = abs(dr/r); end
            if s ~= 0, ea_s = abs(ds/s); end
        else
            % Too few valid equations to solve DR and DS.
            % Need a new guess of R and S to restart the iteration.
            r = r + 1; s = s + 1;
            iter_rs = 0;
        end
    end
    
    fac(iter_fac,2:end) = [-r,-s];  % X^2 + R*X + S*X
    [root(n),root(n-1)] = findQuadRoot(fac(iter_fac,:));
    n = n - 2;
    A = B(3:end);
end

remainder = flip(A);

if n <= 0
    error('Error: Input coefficient vector should include at least 2 elements.')
elseif n == 1  % A(1)*X + A(0) = 0 -> X = -A(0) / A(1)
    root(n) = -A(1)/A(2);
elseif n == 2
    [root(n),root(n-1)] = findQuadRoot([A(3),A(2),A(1)]);
end
