function k=chernoff_bound(b,mu_1,v_1,mu_2,v_2,prior_1)
% Log (base 10) of Chernoff upper bound of classification error between two normals.
% Author:
%   Abhranil Das <abhranil.das@utexas.edu>
%	Center for Perceptual Systems, University of Texas at Austin
% If you use this code, please cite:
%   A new method to compute classification error
%   https://jov.arvojournals.org/article.aspx?articleid=2750251

k=b*log(prior_1) + (1-b)*log(1-prior_1) -...
    b*(1-b)/2*((mu_2-mu_1)'*inv(b*v_1+(1-b)*v_2)*(mu_2-mu_1))...
    + log(det(b*v_1+(1-b)*v_2)/(det(v_1)^b*det(v_2)^(1-b)))/2;

% convert to base 10:
k=k/log(10);
