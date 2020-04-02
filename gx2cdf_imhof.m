function p=gx2cdf_imhof(x,lambda,m,delta,tail)
% Syntax:
% p=gx2cdf_imhof(x,lambda,m,delta)
% p=gx2cdf_imhof(x,lambda,m,delta,'tail')

% Description:
% Returns the CDF of a generalized chi-squared (a weighted sum of
% non-central chi-squares), using Imhof's [1961] algorithm.

% Example:
% [p,pc]=gx2cdf_imhof(25,[1 -5 2],[1 2 3],[2 3 7])

% Inputs:
% x         point at which to evaluate the CDF
% lambda    row vector of coefficients of the non-central chi-squares
% m         row vector of degrees of freedom of the non-central chi-squares
% delta     row vector of non-centrality paramaters (sum of squares of
%           means) of the non-central chi-squares

% Outputs:
% p         computed CDF

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
        f=sin(theta)./(u.*rho);
    end

if nargin==4
    % compute the integral
    p=0.5-integral(@(u) imhof_integrand(u,x,lambda',m',delta'),0,inf)/pi;
    
    if (p<1e-3)||(p>1-1e-3)
        warning("Tail probability may be inaccurate, and might be improved with the 'tail' flag.")
    end
    
elseif nargin>4 % compute tail approximations
    j=(1:3)';
    c=sum((lambda.^j).*(j.*delta+m),2);
    h=c(2)^3/c(3)^2;    
    if c(3)>0
        y=(x-c(1))*sqrt(h/c(2))+h;
        if strcmpi(tail,'lower')
            p=chi2cdf(y,h);
        elseif strcmpi(tail,'upper')
            p=chi2cdf(y,h,'upper');
        end
    else
        c=sum(((-lambda).^j).*(j.*delta+m),2);
        y=(-x-c(1))*sqrt(h/c(2))+h;
        if strcmpi(tail,'lower')
            p=chi2cdf(y,h,'upper');
        elseif strcmpi(tail,'upper')
            p=chi2cdf(y,h);
        end
    end
    
end

end