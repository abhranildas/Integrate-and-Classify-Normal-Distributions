function [init_sign,x]=opt_norm_reg(n,norm_1,norm_2,varargin)
% Return distances and signs of the optimal boundary between normals a and 
% b (relative to mean 1, or optional origin), in the direction
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
addRequired(parser,'norm_1',@(x) isnumeric(x));
addRequired(parser,'norm_2',@(x) isnumeric(x));
addParameter(parser,'orig',norm_1(:,1),@(x) isnumeric(x));
addParameter(parser,'prior_1',0.5, @(x) isnumeric(x) && isscalar(x) && (x > 0) && (x < 1));
addParameter(parser,'vals',eye(2), @(x) isnumeric(x) && ismatrix(x));
parse(parser,n,norm_1,norm_2,varargin{:});

% parse inputs
mu_1=norm_1(:,1);
v_1=norm_1(:,2:end);
mu_2=norm_2(:,1);
v_2=norm_2(:,2:end);
orig=parser.Results.orig;
priors(1)=parser.Results.prior_1;
priors(2)=1-priors(1);
vals=parser.Results.vals;

% optimal quadratic boundary coefficients
[a2,a1,a0]=opt_norm_reg_quad(mu_1,v_1,mu_2,v_2,vals,priors);
coeffs.a2=a2;
coeffs.a1=a1;
coeffs.a0=a0;

% find boundary distance and sign
[init_sign,x]=quad_reg(n,coeffs,orig);