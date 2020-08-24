function [merged_init_sign,merged_x]=opt_class_multi(n,mus,vs,idx,varargin)
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
n_normals=size(mus,2);
parser = inputParser;
addRequired(parser,'n',@isnumeric);
addRequired(parser,'mus',@isnumeric);
addRequired(parser,'vs',@isnumeric);
addRequired(parser,'idx',@isnumeric);
addParameter(parser,'mu',mus(:,idx),@isnumeric);
addParameter(parser,'v',eye(dim),@isnumeric);
addParameter(parser,'priors',ones(1,n_normals)/n_normals, @(x) isnumeric(x) && all(x > 0) && all(x < 1));
addParameter(parser,'vals',eye(n_normals), @(x) isnumeric(x) && ismatrix(x));
parse(parser,n,mus,vs,idx,varargin{:});
%idx=parser.Results.idx;
mu=parser.Results.mu;
v=parser.Results.v;
if isempty(mu) % if origin is not given,
    mu=mus(:,idx); % take the mean of the normal whose boundary is being computed
end
vals=parser.Results.vals;
priors=parser.Results.priors;

other_idxs=[1:idx-1, idx+1:n_normals];
reglist=cell(length(other_idxs),1);
for i=1:length(reglist)
    other_idx=other_idxs(i);
    reglist{i}=@(n,mu,v) ray_scan(opt_class_quad([mus(:,idx), vs(:,:,idx)],[mus(:,other_idx), vs(:,:,other_idx)],...
        'prior_1',priors(idx)/(priors(idx)+priors(other_idx)),'vals',vals([idx other_idx],[idx other_idx])),n,'mu',mu,'v',v);
end

[merged_init_sign,merged_x]=combine_regs(reglist,'and',n,'mu',mu,'v',v);
end