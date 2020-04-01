function Q=gx2cdf_imhof(x,lambda,m,delta)
% Returns the CDF of a generalized chi-squared (a weighted sum of
% non-central chi-squares), using Imhof's [1961] algorithm.

% Inputs:
% x         point at which to evaluate the CDF
% lambda    row vector of coefficients of the non-central chi-squares
% m         row vector of degrees of freedom of the non-central chi-squares
% delta     row vector of non-centrality paramaters (sum of squares of
%           means) of the non-central chi-squares

% Outputs:
% Q         computed CDF

% Example:
% Q=gx2cdf_imhof(25,[1 -5 2],[1 2 3],[2 3 7])

% Author:
% Abhranil Das <abhranil.das@utexas.edu>
% Center for Perceptual Systems, University of Texas at Austin

% If you use this code, you may cite:
% A new method to compute classification error
% jov.arvojournals.org/article.aspx?articleid=2750251

    % define the integrand (lambda, m, delta must be column vectors here)
    function f=imhof_integrand(u,x,lambda,m,delta)
        theta=sum(m.*atan(lambda*u)+(delta.*(lambda*u))./(1+lambda.^2*u.^2),1)/2-u*x/2;
        rho=prod(((1+lambda.^2*u.^2).^(m/4)).*exp(((lambda.^2*u.^2).*delta)./(2*(1+lambda.^2*u.^2))),1);
        f=exp(-u)/2-sin(theta)./(pi*u.*rho);
    end

% compute the integral
Q=integral(@(u) imhof_integrand(u,x,lambda',m',delta'),0,inf);

if Q<0 % integrals smaller than precision limit may return negative
    Q=nan;
end
end