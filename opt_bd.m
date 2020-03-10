function [r,r_sign]=opt_bd(n,dist_a,dist_b,varargin)
% Return distances and signs of the optimal boundary between normals a and 
% b (standardized wrt normal a, or optional normal_wrt), in the direction
% of vector(s) n.
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

% parse inputs
parser = inputParser;
addRequired(parser,'n',@(x) isnumeric(x));
addRequired(parser,'dist_a',@(x) isnumeric(x));
addRequired(parser,'dist_b',@(x) isnumeric(x));
addParameter(parser,'dist_wrt',[],@(x) isnumeric(x));
addParameter(parser,'prior_a',0.5, @(x) isnumeric(x) && isscalar(x) && (x > 0) && (x < 1));
addParameter(parser,'vals',eye(2), @(x) isnumeric(x) && ismatrix(x));
parse(parser,n,dist_a,dist_b,varargin{:});

% parse inputs
mu_a=dist_a(:,1);
v_a=dist_a(:,2:end);
mu_b=dist_b(:,1);
v_b=dist_b(:,2:end);
dist_wrt=parser.Results.dist_wrt;
if isempty(dist_wrt)
    mu_wrt=mu_a;
    v_wrt=v_a;
else
    mu_wrt=dist_wrt(:,1);
    v_wrt=dist_wrt(:,2:end);
end
prior_a=parser.Results.prior_a;
prior_b=1-prior_a;
vals=parser.Results.vals;

% optimal quadratic boundary coefficients
coeffs.a2=inv(v_b)-inv(v_a);
coeffs.a1=2*(v_a\mu_a-v_b\mu_b);
coeffs.a0=mu_b'/v_b*mu_b-mu_a'/v_a*mu_a+log((((vals(1,1)-vals(1,2))*prior_a)/((vals(2,2)-vals(2,1))*prior_b))^2*det(v_b)/det(v_a));

% find boundary distance and sign
[r,r_sign]=quad_bd(n,mu_wrt,v_wrt,coeffs);