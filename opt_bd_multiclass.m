function [r,r_sign]=opt_bd_multiclass(n,normals,varargin)
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
parser = inputParser;
addParameter(parser,'vals',[], @(x) isnumeric(x) && ismatrix(x)); % TODO default vals
parse(parser,varargin{:});


mu_a=dist_a(:,1);
v_a=dist_a(:,2:end);
mu_b=dist_b(:,1);
v_b=dist_b(:,2:end);

p_a=parser.Results.p_a;
p_b=1-p_a;
vals=parser.Results.vals;

% find boundary distance and sign
[r,r_sign]=quad_bd(n,mu_a,v_a,coeffs);