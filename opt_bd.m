function [r,r_sign]=opt_bd(n,dist_a,dist_b,varargin)
% Return distances and signs (relative to mu_a) of the optimal quadratic
% boundary between normals a and b, in the unit direction n.
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
p = inputParser;
addRequired(p,'n',@(x) isnumeric(x));
addRequired(p,'dist_a',@(x) isnumeric(x));
addRequired(p,'dist_b',@(x) isnumeric(x));
addParameter(p,'p_a',0.5, @(x) isnumeric(x) && isscalar(x) && (x > 0) && (x < 1));
addParameter(p,'vals',eye(2), @(x) isnumeric(x) && ismatrix(x));
parse(p,n,dist_a,dist_b,varargin{:});

mu_a=dist_a(:,1);
v_a=dist_a(:,2:end);
mu_b=dist_b(:,1);
v_b=dist_b(:,2:end);

p_a=p.Results.p_a;
p_b=1-p_a;
vals=p.Results.vals;

% optimal quadratic boundary coefficients
coeffs.a2=inv(v_b)-inv(v_a);
coeffs.a1=2*(v_a\mu_a-v_b\mu_b);
coeffs.a0=mu_b'/v_b*mu_b-mu_a'/v_a*mu_a+log((((vals(1,1)-vals(1,2))*p_a)/((vals(2,2)-vals(2,1))*p_b))^2*det(v_b)/det(v_a));

% find boundary distance and sign
[r,r_sign]=quad_bd(n,mu_a,v_a,coeffs);