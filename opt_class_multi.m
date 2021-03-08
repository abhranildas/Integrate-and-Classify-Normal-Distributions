function [merged_init_sign,merged_x,merged_samp_correct]=opt_class_multi(n,normals,idx,varargin)
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
n_normals=length(normals);
parser=inputParser;
parser.KeepUnmatched=true;
addRequired(parser,'n',@isnumeric);
addRequired(parser,'normals',@isstruct);
addRequired(parser,'idx',@isnumeric);
addParameter(parser,'orig',normals(idx).mu,@isnumeric);
addParameter(parser,'priors',ones(1,n_normals)/n_normals, @(x) isnumeric(x) && all(x > 0) && all(x < 1));
addParameter(parser,'vals',eye(n_normals), @(x) isnumeric(x) && ismatrix(x));
parse(parser,n,normals,idx,varargin{:});
orig=parser.Results.orig;
vals=parser.Results.vals;
priors=parser.Results.priors;

other_idxs=[1:idx-1, idx+1:n_normals];
domlist=cell(length(other_idxs),1);
for i=1:length(domlist)
    other_idx=other_idxs(i);
    domlist{i}=@(n,orig) standard_ray_scan(opt_class_quad([normals(idx).mu, normals(idx).v],[normals(other_idx).mu, normals(other_idx).v],...
        'prior_1',priors(idx)/(priors(idx)+priors(other_idx)),'vals',vals([idx other_idx],[idx other_idx])),n,'mu',orig);
end

if nargout==3
    [~,~,merged_samp_correct]=combine_ray_scan_doms(domlist,'and',n,'orig',orig);
    merged_init_sign=[];
    merged_x=[];
else
    [merged_init_sign,merged_x]=combine_ray_scan_doms(domlist,'and',n,'orig',orig);
end
end