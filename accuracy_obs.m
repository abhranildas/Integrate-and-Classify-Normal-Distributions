function [acc_obs, acc_obs_a, acc_obs_b]=accuracy_obs(obs_a,obs_b,bd_coeffs)

a2=bd_coeffs.a2;
a1=bd_coeffs.a1;
a0=bd_coeffs.a0;

correct_obs_a=dot(obs_a,obs_a*a2',2) + obs_a*a1 + a0 > 0;
correct_obs_b=dot(obs_b,obs_b*a2',2) + obs_b*a1 + a0 < 0;

acc_obs_a=mean(correct_obs_a);
acc_obs_b=mean(correct_obs_b);
acc_obs=mean([correct_obs_a;correct_obs_b]);

end


