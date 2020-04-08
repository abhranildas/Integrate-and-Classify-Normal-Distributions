function [init_sign,x,samp_correct]=invert_reg(reg,n,orig)
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

if nargout==3
    [~,~,samp_correct]=reg(n,orig);
    samp_correct=~samp_correct;
    init_sign=[];
    x=[];    
else    
    [init_sign,x]=reg(n,orig);
    init_sign=-init_sign;
end