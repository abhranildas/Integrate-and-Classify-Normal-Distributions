%% 2D example
mu_a=[0;0];
sig_a=[[1 0]; 
       [0 1]];

mu_b=[0;0];
sig_b=[[2 0]; 
       [0 2]];
p_a=.5;

[dPrime,err_rate]=d_prime(mu_a,sig_a,mu_b,sig_b,p_a,1)

%% 3D example
mu_a=[0;0;0];
sig_a=[[1 0 0]; 
       [0 1 0];
       [0 0 1]];

mu_b=[0;0;1];
sig_b=[[1 0 0]; 
       [0 1 0];
       [0 0 1]];
p_a=.5;

[dPrime,err_rate]=d_prime(mu_a,sig_a,mu_b,sig_b,p_a,1)