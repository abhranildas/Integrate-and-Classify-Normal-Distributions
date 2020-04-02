function [p,pc]=int_stdnorm_quad_gx2(coeffs)
% Find the probability that a quadratic form of a standard normal variate z
% z'a2z + a1'z + a0 >= 0
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

a2=coeffs.a2;
a1=coeffs.a1;
a0=coeffs.a0;

if ~nnz(a2) % if a2 is zero, linear discriminant
    p=normcdf(a0/norm(a1));
    pc=normcdf(-a0/norm(a1)); % complement of p. It's useful to return it when small, and p is rounded to 1.
else
    [R,D]=eig(a2);
    d=diag(D);
    b2=(R*a1).^2; % b^2
    c=a0-sum(b2./(4*d));
    
    [lambda,~,ic]=unique(d); % unique eigenvalues
    m=accumarray(ic,1); % total dof of each eigenvalue
    delta=arrayfun(@(x) sum(b2(d==x)),lambda)./(4*lambda.^2); % total non-centrality for each eigenvalue
    
    % use Imhof's method to compute the CDF.
    p_lower=gx2cdf_imhof(-c,lambda',m',delta','lower');
    p_upper=gx2cdf_imhof(-c,lambda',m',delta','upper');
    pc_actual=gx2cdf_imhof(-c,lambda',m',delta');
    if min(p_lower,p_upper)>1e-3 % use actual Imhof's method
        pc=pc_actual;
        p=1-pc;
    else % use tail approximation
        warning('Using the tail approximation of the generalized chi-squared CDF.')
        if p_lower<p_upper % lower tail
            pc=p_lower;
            p=1-pc;
        else % upper tail
            p=p_upper;
            pc=1-p;
        end
    end
end
