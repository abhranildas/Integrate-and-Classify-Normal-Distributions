function [merged_init_sign,merged_x]=combine_regs(reglist,op,n,orig)
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

n_regs=length(reglist);
n_dirs=size(n,2);

% compile distances and signs to all boundaries
all_init_sign=nan(n_regs,n_dirs); all_x=cell(1,n_dirs);
for i=1:n_regs
    [init_sign,x]=reglist{i}(n,orig);
    all_x=arrayfun(@(x,y) [x{:};y], all_x,x,'un',0);
    all_init_sign(i,:)=init_sign;
end

[merged_init_sign,merged_x]=cellfun(@(all_init_sign_ray,all_x_ray) combine_regs_ray(all_init_sign_ray,all_x_ray,op),num2cell(all_init_sign,1),all_x,'un',0);
merged_init_sign=cell2mat(merged_init_sign);
end