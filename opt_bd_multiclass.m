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

%vals=parser.Results.vals;

% find boundary distance and sign
[r,r_sign]=opt_bd(n,[normals(1).mu, normals(1).v],[normals(2).mu, normals(2).v],...
    'p_a',normals(1).p/(normals(1).p+normals(2).p));