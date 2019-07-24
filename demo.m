% Examples for classify library.
% Author: Abhranil Das <abhranil.das@utexas.edu>
% Please cite if you use this code.

%% 1D, unequal prior
mu_a=0;
v_a=1;

mu_b=0.5;
v_b=1.5;

results=bayes_classify([mu_a,v_a],[mu_b,v_b],'p_a',.7)

%% 1D, input observations, unequal occurrences
mu_a=0;
v_a=1;
obs_a=normrnd(mu_a,sqrt(v_a),[700 1]);

mu_b=0.5;
v_b=1.5;
obs_b=normrnd(mu_b,sqrt(v_b),[300 1]);

results=bayes_classify(obs_a,obs_b,'type','obs')

%% 2D, simple

mu_a=[4; 5];
v_a=[2 1.1; 1.1 1];

mu_b=[4; 5];
v_b=1.4*v_a;

results=bayes_classify([mu_a,v_a],[mu_b,v_b]);

%% 2D, input observations

mu_a=[2 4];
v_a=[1 1.5; 1.5 3];
obs_a = mvnrnd(mu_a,v_a,1000);

mu_b=[5 0];
v_b=[3 0; 0 1];
obs_b = mvnrnd(mu_b,v_b,1000);

results1=bayes_classify(obs_a,obs_b,'type','obs');

axis image
xlim([-10 10])
ylim([-10 10])

% modify boundary
custom_bd_coeffs=results1.bd_coeffs_obs_opt;
custom_bd_coeffs.a2=custom_bd_coeffs.a2+.2;
custom_bd_coeffs.a1=custom_bd_coeffs.a1-5;
custom_bd_coeffs.a0=custom_bd_coeffs.a0+10;

results2=bayes_classify(obs_a,obs_b,'type','obs','custom_bd_coeffs',custom_bd_coeffs);

axis image
xlim([-10 10])
ylim([-10 10])

%% 2D, input non-Gaussian observations
n_samp=100;
obs_a=unifrnd(0,3,[n_samp 2]);
obs_b=unifrnd(1,4,[n_samp 2]);

results=bayes_classify(obs_a,obs_b,'type','obs')
xlim([-1 5])
ylim([-1 5])

%% 3D, simple

mu_a=[0;0;0];
v_a=eye(3);

mu_b=[2;1;1];
v_b=2*eye(3);

results=bayes_classify([mu_a,v_a],[mu_b,v_b]);

%% 3D, far apart
mu_a=[0;0;0];
v_a=eye(3);

mu_b=[100;0;0];
v_b=1.5*eye(3);

results=bayes_classify([mu_a,v_a],[mu_b,v_b])

%% 3D, complex

mu_a=[0;0;0];
v_a=[1 .5 0; .5 2 1; 0 1 4];

mu_b=[2; 1; 2];
v_b=[2 1 -1.5; 1 3 -2; -1.5 -2 4];

results=bayes_classify([mu_a,v_a],[mu_b,v_b]);

%% 4D, simple

mu_a=[0;0;0;0];
v_a=eye(4);

mu_b=[1;1;1;1];
v_b=eye(4);

results=bayes_classify([mu_a,v_a],[mu_b,v_b]);

%% 4D, input observations

mu_a=[0;0;0;0];
v_a=eye(4);
obs_a = mvnrnd(mu_a,v_a,1000);

mu_b=[1;1;1;1];
v_b=eye(4);
obs_b = mvnrnd(mu_b,v_b,1000);

results=bayes_classify(obs_a,obs_b,'type','obs');

