function [p,pc]=int_norm_quad_gx2(mu,v,quad,varargin)
% Find the probability that a quadratic form of a normal variate x
% x'q2x + q1'x + q0 >= 0
% using the generalized chi-squared CDF (Imhof's method).
%
% How to use this command:
% See github readme at https://github.com/abhranildas/classify
%
% Credits:
%   Abhranil Das <abhranil.das@utexas.edu>
%	Wilson S Geisler
%	Center for Perceptual Systems, University of Texas at Austin
% If you use this code, please cite:
%   A new method to compute classification error
%   https://jov.arvojournals.org/article.aspx?articleid=2750251

parser = inputParser;
addRequired(parser,'mu',@isnumeric);
addRequired(parser,'v',@isnumeric);
addRequired(parser,'quad');
addParameter(parser,'AbsTol',1e-10);
addParameter(parser,'RelTol',1e-2);
parse(parser,mu,v,quad,varargin{:});

AbsTol=parser.Results.AbsTol;
RelTol=parser.Results.RelTol;

if ~nnz(quad.q2) % if q2 is zero, linear discriminant
    
    % standardize quadratic coefficients
    quad_s=standard_quad(quad,mu,v); % standardize quad
%     q1=sqrtm(v)*(2*quad.q2*mu+quad.q1);
%     q0=mu'*quad.q2*mu+quad.q1'*mu+quad.q0;

    p=normcdf(quad_s.q0/norm(quad_s.q1));
    pc=normcdf(-quad_s.q0/norm(quad_s.q1));
else
    
    % get generalized chi-squared parameters
    [lambda,m,delta,sigma,c]=gx2_params_norm_quad(mu,v,quad);
    
    if (AbsTol==1e-10)&&(RelTol==1e-2)
        pc=gx2cdf(0,lambda,m,delta,sigma,c);
        p=gx2cdf(0,lambda,m,delta,sigma,c,'upper');
    else
        pc=gx2cdf(0,lambda,m,delta,sigma,c,'AbsTol',AbsTol,'RelTol',RelTol);
        p=gx2cdf(0,lambda,m,delta,sigma,c,'upper','AbsTol',AbsTol,'RelTol',RelTol);
    end
end
