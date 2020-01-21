function k=chernoff_bound(b,mu_a,v_a,mu_b,v_b,p_a)
% Log of Chernoff upper bound of classification error between two Gaussians.
% Credits:
%   Abhranil Das <abhranil.das@utexas.edu>
%	R Calen Walshe
%	Wilson S Geisler
%	Center for Perceptual Systems, University of Texas at Austin
% If you use this code, please cite:
%   A new method to compute classification error
%   https://jov.arvojournals.org/article.aspx?articleid=2750251

% k=p_a^b*(1-p_a)^(1-b)* exp(-...
%     [b*(1-b)/2*((mu_b-mu_a)'*inv(b*v_a+(1-b)*v_b)*(mu_b-mu_a))...
%     + log(det(b*v_a+(1-b)*v_b)/(det(v_a)^b*det(v_b)^(1-b)))/2] ...
%     );

k=b*log(p_a) + (1-b)*log(1-p_a) -...
    b*(1-b)/2*((mu_b-mu_a)'*inv(b*v_a+(1-b)*v_b)*(mu_b-mu_a))...
    + log(det(b*v_a+(1-b)*v_b)/(det(v_a)^b*det(v_b)^(1-b)))/2 ...
    ;

% convert to base 10:
k=k/log(10);
