function [merged_init_sign,merged_x]=opt_reg_multi(n,mus,vs,varargin)
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
n_normals=size(mus,2);
parser = inputParser;
addRequired(parser,'n',@(x) isnumeric(x));
addRequired(parser,'mus',@(x) isnumeric(x));
addRequired(parser,'vs',@(x) isnumeric(x));
addParameter(parser,'idx',1, @(x) isnumeric(x));
addParameter(parser,'orig',[],@(x) isnumeric(x));
addParameter(parser,'priors',ones(1,n_normals)/n_normals, @(x) isnumeric(x) && all(x > 0) && all(x < 1));
addParameter(parser,'vals',eye(n_normals), @(x) isnumeric(x) && ismatrix(x));
parse(parser,n,mus,vs,varargin{:});
idx=parser.Results.idx;
orig=parser.Results.orig;
if isempty(orig) % if origin is not given,
    orig=mus(:,idx); % take the mean of the normal whose boundary is being computed
end
vals=parser.Results.vals;
priors=parser.Results.priors;

other_idxs=[1:idx-1, idx+1:n_normals];
reglist=cell(length(other_idxs),1);
for i=1:length(reglist)
    other_idx=other_idxs(i);
    reglist{i}=@(n,orig) ray_scan(opt_reg_quad([mus(:,idx), vs(:,:,idx)],[mus(:,other_idx), vs(:,:,other_idx)],...
        'prior_1',priors(idx)/(priors(idx)+priors(other_idx)),'vals',vals([idx other_idx],[idx other_idx])),'quad',n,orig);
end

[merged_init_sign,merged_x]=combine_regs(reglist,'and',n,orig);
end