%% 2D example
mu_a=[2;3];
sig_a=[[.5 -.2]; 
       [-.2 1]];

mu_b=[4;4];
sig_b=[[3 1]; 
       [1 2]];
p_a=.2;

[dPrime,err_rate]=class_error(mu_a,sig_a,mu_b,sig_b,p_a,1)

%% 3D example
mu_a=[0;0;0];
sig_a=[[1 0 0]; 
       [0 1 0];
       [0 0 1]];

mu_b=[1;0;0];
sig_b=[[1 0 0]; 
       [0 1 0];
       [0 0 1]];
p_a=.5;

custom_dec_coeffs=[];

[dPrime,err_rate,dec_bd_coeffs]=class_error(mu_a,sig_a,mu_b,sig_b,p_a,custom_dec_coeffs,1)