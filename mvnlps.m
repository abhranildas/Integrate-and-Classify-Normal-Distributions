function mvnval  = mvnlps( mu, sigma, q, e, r, re )
%
% MVNLPS Multivariate Normal Distribution Value for an ellipsoid.        
%        MVNVAL = MVNLPS( MU, SIGMA, Q, E, R, RE ) computes
%        an MVN value to relative accuracy RE for an ellipsoid centered 
%        at Q with radius R and positive semi-definite ellipsoid matrix E: 
%                MVNVAL = PROB( ( X - Q )'E ( X - Q ) < R^2 )
%        SIGMA is a positive definite covariance matrix for a multivariate
%        normal (MVN) density with mean MU. MU and Q must be column vectors.
%        Example:
%         sg = [ 3 2 1;2 2 1;1 1 1]; mu = [1 -1 0]';        
%         e = [4 1 -1; 1 2 1; -1 1 2]; q = [2 3 -2]';
%         p = mvnlps( mu, sg, q, e, 4, 1e-5 ); disp(p)
%
%        Reference for basic algorithm is (S&O):
%         Sheil, J. and O'Muircheartaigh, I. (1977), Algorithm AS 106:
%           The Distribution of Non-negative Quadratic Forms in Normal
%           Variables, Applied Statistics 26, pp. 92-98.
%
%        Matlab implementation by
%         Alan Genz, WSU Math, Pullman, WA 99164-3113; alangenz@wsu.edu.
%
%
% Copyright (C) 2013, Alan Genz,  All rights reserved.               
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided the following conditions are met:
%   1. Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%   2. Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in 
%      the documentation and/or other materials provided with the 
%      distribution.
%   3. The contributor name(s) may not be used to endorse or promote 
%      products derived from this software without specific prior 
%      written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
% FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
% COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
% INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
% BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS 
% OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND 
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR 
% TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%  Transformation to problem with diagonal covariance matrix
%
   [ n, n ] = size(sigma); L = chol(sigma); [ V D ]  = eig( L*e*L' );
   cov = diag(D);  d = sqrt(D)*( inv(L)*V )'*( q - mu );
%
%  Basic S&O algorithm follows
%
   kmx = 2000; covmx = max(cov); rsqrd = r*r; np = 0;
   for j = 1 : n
     if cov(j) > 1e-10*covmx, np = np + 1; 
       alpha(np) = cov(j); bsqrd(np) = d(j)^2/alpha(np);       
     else, rsqrd = rsqrd - d(j)^2; 
     end
   end, A = 1/prod(alpha); lambda = sum(bsqrd); gam = ones(1,np); 
   if rsqrd <= 0, mvnval = 0; 
   else, bet = 29*min(alpha)/32; tbeta = rsqrd/bet; 
     A = sqrt(A*bet^np); betalph = bet./alpha;
%
%     c(1), c(2), ..., are S&O's c_0, c_1, ...;
%     bsqrd(j) is S&O's b_j^2 ; tbeta is S&O's t/beta, np is S&O's n';
%     g(j) is S&O's g_j; gam(j) is S&O's gamma_j^k.  
%
      csum = A*exp(-lambda/2); c(1) = csum; lgb = log(tbeta);  
      if mod( np, 2 )  == 1 
         F = erf(sqrt(tbeta/2)); lgf = -( tbeta + log(2*pi*tbeta) )/2;
         for nc = 3:2:np, lgf = lgb + lgf - log(nc-2); F = F - 2*exp(lgf); end
      else  
         F = 1 - exp(-tbeta/2);  lgf = - tbeta/2 - log(2); 
         for nc = 4:2:np, lgf = lgb + lgf - log(nc-2); F = F - 2*exp(lgf); end
      end
      % equivalent F = gammainc( tbeta/2 , np/2 );  
      mvnval = c(1)*F; 
%
%     for-loop computes S&O series  
%
      for k = 1 : kmx 
	g(k) = gam*( 1+betalph.*(k*bsqrd-1))'/2; gam = gam.*(1-betalph);
	c(k+1) = g(k:-1:1)*c(1:k)'/k; csum = csum + c(k+1);
	lgf = lgb + lgf - log(np+2*k-2); F = F - 2*exp(lgf);
	% equivalent F = gammainc( tbeta/2, np/2 + k ); 
	mvnval = mvnval + c(k+1)*F; 
	if ( 1 - csum )*F < re*mvnval, break; end
      end
   end
%
% End function mvnlps
%
