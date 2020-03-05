function [ex_val_obs, ex_val_obs_a, ex_val_obs_b, outcome_vals_obs]=val_obs(obs_a,obs_b,bd_coeffs,vals)
% Expected value given two samples and a boundary. If outcome values
% are not additionally specified, this is the classification accuracy.
% Credits:
%   Abhranil Das <abhranil.das@utexas.edu>
%	R Calen Walshe
%	Wilson S Geisler
%	Center for Perceptual Systems, University of Texas at Austin
% If you use this code, please cite:
%   A new method to compute classification error
%   https://jov.arvojournals.org/article.aspx?articleid=2750251

a2=bd_coeffs.a2;
a1=bd_coeffs.a1;
a0=bd_coeffs.a0;

correct_obs_a=dot(obs_a,obs_a*a2',2) + obs_a*a1 + a0 > 0;
correct_obs_b=dot(obs_b,obs_b*a2',2) + obs_b*a1 + a0 < 0;

outcome_vals_obs=zeros(2);
outcome_vals_obs(1,1)=sum(correct_obs_a)*vals(1,1);
outcome_vals_obs(1,2)=sum(~correct_obs_a)*vals(1,2);
outcome_vals_obs(2,2)=sum(correct_obs_b)*vals(2,2);
outcome_vals_obs(2,1)=sum(~correct_obs_b)*vals(2,1);

ex_val_obs_a=sum(outcome_vals_obs(1,:))/size(obs_a,1);
ex_val_obs_b=sum(outcome_vals_obs(2,:))/size(obs_b,1);

ex_val_obs=sum(outcome_vals_obs(:))/(size(obs_a,1)+size(obs_b,1));

% acc_obs_a=mean(correct_obs_a);
% acc_obs_b=mean(correct_obs_b);
% acc_obs=mean([correct_obs_a;correct_obs_b]);

end


