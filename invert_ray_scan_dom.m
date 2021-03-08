function [init_sign,x,samp_correct]=invert_ray_scan_dom(dom,n,varargin)
% Return distances and signs of the optimal boundary between normal 1 and
% several others (standardized wrt normal 1, or optional normal_wrt), in the direction
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
dim=size(n,1);
parser = inputParser;
addRequired(parser,'reg');
addRequired(parser,'n',@isnumeric);
addParameter(parser,'orig',zeros(dim,1),@isnumeric);

parse(parser,dom,n,varargin{:});
orig=parser.Results.orig;

if nargout==3
    [~,~,samp_correct]=dom(n,[]);
    samp_correct=~samp_correct;
    init_sign=[];
    x=[];    
else    
    [init_sign,x]=dom(n,orig);
    init_sign=-init_sign;
end