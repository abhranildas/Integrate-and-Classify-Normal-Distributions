function [lambda,m,delta,c]=norm_quad_to_gx2_params(mu,v,quad)
% A quadratic form of a normal variable is distributed as a generalized
% chi-squared. This function takes the normal parameters and the quadratic
% coeffs and returns the parameters of the generalized chi-squared.

% standardize the space
a2=sqrtm(v)*quad.a2*sqrtm(v);
a1=sqrtm(v)*(2*quad.a2*mu+quad.a1);
a0=mu'*quad.a2*mu+quad.a1'*mu+quad.a0;

[R,D]=eig(a2);
d=diag(D)';
b2=(a1'*R').^2; % b^2

[lambda,~,ic]=unique(d); % unique eigenvalues
m=accumarray(ic,1)'; % total dof of each eigenvalue
delta=arrayfun(@(x) sum(b2(d==x)),lambda)./(4*lambda.^2); % total non-centrality for each eigenvalue
c=a0-sum(b2./(4*d));
