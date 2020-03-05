function [r,r_sign]=quad_bd(n,mu,v,coeffs)
% Return distances and signs of a standardized quadratic boundary in
% the direction(s) of vector(s) n. For n=0, returns if mu is in the region.
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

C=chol(v,'lower'); % Cholesky decomposition of vcov

% standardize boundary coefficients
a2=C'*coeffs.a2*C;
a1=C'*(2*coeffs.a2*mu+coeffs.a1);
a0=mu'*coeffs.a2*mu+coeffs.a1'*mu+coeffs.a0;

if all(~n) % if at the origin mu
    r=nan;
    r_sign=double(a0>0); % return whether mu is in positive region
else
    % normalize direction vector
    n=n./vecnorm(n);
    
    % find boundary distance and sign
    %q2=n'*a2*n;
    q2=dot(n,a2*n);
    q1=a1'*n;
    %r=roots([q2 q1 a0])';
    r=arrayfun(@(x,y) roots([x y a0])',q2,q1,'un',0);
    r=cellfun(@(x) x(~imag(x)),r,'un',0); % only real roots
    %r=r(~imag(r)); % only real roots
    %rs=rs(rs>0); % only positive roots
    r_sign=cellfun(@(x,a,b) sign(2*a*x.^2+b*x),...
        r,num2cell(q2),num2cell(q1),'un',0);
    %r_list=horzcat(r_list{:});    
    %r_sign_list=horzcat(r_sign_list{:});
    
    % only real roots:
    %r_real_list=r_list(~imag(r_list));
    %r_real_sign_list=r_sign_list(~imag(r_list));
    %r_sign=sign(2*q2.*r.^2+q1.*r);
end