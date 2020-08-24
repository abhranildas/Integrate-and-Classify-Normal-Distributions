function [lambda,m,delta,sigma,c]=gx2_params_norm_quad(mu,v,quad)
% A quadratic form of a normal variable is distributed as a generalized
% chi-squared. This function takes the normal parameters and the quadratic
% coeffs and returns the parameters of the generalized chi-squared.

% Example:

% mu=[1;2]; % mean
% v=[2 1; 1 3]; % covariance matrix

% Say q(x)=(x1+x2)^2-x1-1 = [x1;x2]'*[1 1; 1 1]*[x1;x2] + [-1;0]'*[x1;x2] - 1:
% quad.q2=[1 1; 1 1];
% quad.q1=[-1;0];
% quad.q0=-1;

% [lambda,m,delta,sigma,c]=gx2_params_norm_quad(mu,v,quad)

% Inputs:
% mu        column vector of normal mean
% v         normal covariance matrix
% quad      struct with quadratic form coefficients:
%               q2      matrix of quadratic coefficients
%               q1      column vector of linear coefficients
%               q0      scalar constant

% Output:
% lambda    row vector of coefficients of the non-central chi-squares
% m         row vector of degrees of freedom of the non-central chi-squares
% delta     row vector of non-centrality paramaters (sum of squares of
%           means) of the non-central chi-squares
% sigma     sd of normal term
% c         constant term

% Author:
% Abhranil Das <abhranil.das@utexas.edu>
% Center for Perceptual Systems, University of Texas at Austin

% If you use this code, you may cite:
% A new method to compute classification error
% jov.arvojournals.org/article.aspx?articleid=2750251


% standardize the space
q2=sqrtm(v)*quad.q2*sqrtm(v);
q1=sqrtm(v)*(2*quad.q2*mu+quad.q1);
q0=mu'*quad.q2*mu+quad.q1'*mu+quad.q0;

[R,D]=eig(q2);
d=diag(D)';
b=(R'*q1)';

[lambda,~,ic]=unique(nonzeros(d)'); % unique non-zero eigenvalues
m=accumarray(ic,1)'; % total dof of each eigenvalue
delta=arrayfun(@(x) sum((b(d==x)).^2),lambda)./(4*lambda.^2); % total non-centrality for each eigenvalue
sigma=norm(b(~d));
c=q0-sum(lambda.*delta);
