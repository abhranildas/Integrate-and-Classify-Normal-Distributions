function acc_obs=accuracy_obs_optimize(dim,x,obs_a,obs_b)
% Function to optimize classification boundary between two samples.
% Author: Abhranil Das <abhranil.das@utexas.edu>
% Please cite if you use this code.

a2=reshape(x(1:dim^2),[dim dim])';
a1=x(dim^2+1:dim^2+dim);
a0=x(end);

acc_obs=-sum([dot(obs_a,obs_a*a2',2) + obs_a*a1 + a0 > 0;...
    dot(obs_b,obs_b*a2',2) + obs_b*a1 + a0 < 0]);