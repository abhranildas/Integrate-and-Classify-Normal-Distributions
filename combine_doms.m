function [merged_init_sign,merged_x,merged_samp_correct]=combine_doms(domlist,op,n,varargin)
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
parser=inputParser;
addRequired(parser,'domlist');
addRequired(parser,'op');
addRequired(parser,'n',@isnumeric);
addParameter(parser,'mu',zeros(dim,1),@isnumeric);
addParameter(parser,'v',eye(dim),@isnumeric);

parse(parser,domlist,op,n,varargin{:});
mu=parser.Results.mu;
v=parser.Results.v;

n_doms=length(domlist);
n_dirs=size(n,2);

if nargout==3
    all_samp_correct=nan(n_doms,n_dirs);
    for i=1:n_doms
        [~,~,samp_correct]=domlist{i}(n,[],[]);
        all_samp_correct(i,:)=samp_correct;
    end
    if strcmpi(op,'or')
        merged_samp_correct=any(all_samp_correct);
    elseif strcmpi(op,'and')
        merged_samp_correct=all(all_samp_correct);
    end
    merged_init_sign=[];
    merged_x=[];
    
else    
    % compile distances and signs to all boundaries
    all_init_sign=nan(n_doms,n_dirs); all_x=cell(1,n_dirs);
    for i=1:n_doms
        [init_sign,x]=domlist{i}(n,mu,v);
        all_x=arrayfun(@(x,y) [x{:};y], all_x,x,'un',0);
        all_init_sign(i,:)=init_sign;
    end
    
    [merged_init_sign,merged_x]=cellfun(@(all_init_sign_ray,all_x_ray) combine_doms_ray(all_init_sign_ray,all_x_ray,op),num2cell(all_init_sign,1),all_x,'un',0);
    merged_init_sign=cell2mat(merged_init_sign);
end