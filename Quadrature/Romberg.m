function [tab,etab] = Romberg(x,f,nrow0,es,nintv0)
%==========================================================================
%ROMBERG return the Romberg quadrature table given the xdata and function.
%
% Author: Ma Chenglong (github: Masaaki-75)
%         Copyright 2021 Sun Yat-sen University.
%
% Syntax:
%
%   [TAB,ETAB] = ROMBERG(X,F,NROW0,ES,NINTV0);
%       
%       X is a vector containing points in the interval of the integral; F
%       is a function handle for integrand; 
%
%       NROW0 is a positive integer specifying the number of rows in the 
%       output Romberg quadrature table; ES represents relative error
%       tolerance; NINTV0 is the number of subintervals.
%
%       The output Romberg table TAB takes the following form:
%
%       h=h0         -> [I[1,1], I[1,2], ..., I[m,n]]
%       h=h0/2       -> [I[2,1], I[2,2], ...,    0  ]
%       h=h0/4       -> [I[3,1], I[3,2], ...,    0  ]
%                       [  :   ,    :  , ...,    0  ]
%       h=h0/2^(m-1) -> [I[m,1],    0  , ...,    0  ]
%
%==========================================================================

if ~exist('es','var')||isempty(es), es = 1e-6; end
if ~exist('nintv0','var')||isempty(nintv0), nintv0 = 1; end
if ~exist('nrow0','var')||isempty(nrow0), nrow0 = 2; end

tab = zeros(nrow0);
etab = ones(nrow0);

for jj = 1:nrow0
    % nintv0*ii is the number of subintervals
    alpha = 2^(jj-1);
    tab(jj,1) = TrapzSubintvN(x,f,nintv0*alpha);
end

for kk = 2:nrow0
    for jj = 1:nrow0-kk+1
        t = 4^(kk-1);
        tab(jj,kk) = t * tab(jj+1,kk-1) - tab(jj,kk-1);
        tab(jj,kk) = tab(jj,kk)/(t-1);
        
        if tab(jj,kk) ~= 0
            etab(jj,kk) = abs(1 - tab(jj+1,kk-1)/tab(jj,kk));
        end
        
        if etab(jj,kk) <= es
            tab = tab(:,1:kk);
            etab = etab(:,1:kk);
            outstr = 'Iteration stops because the result has achieved the required accuracy';
            errstr = [' (',num2str(etab(jj,kk)*100),'% / ',num2str(es*100),'%)'];
            disp([newline,outstr,errstr,newline])
            return;
        end
    end
end

outstr = 'Iteration stops because maximum number of iterations allowed has been reached.';
errstr = [' (',num2str(etab(jj,kk)*100),'% / ',num2str(es*100),'%)'];
disp([newline,outstr,errstr,newline])
