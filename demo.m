%% 1D, unequal prior
mu_a=0;
v_a=1;

mu_b=0.5;
v_b=1.5;

results=classify([mu_a,v_a],[mu_b,v_b],'p_a',.7)

%% 1D, far apart
mu_a=0;
v_a=1;

mu_b=100;
v_b=1.5;

results=classify([mu_a,v_a],[mu_b,v_b])

%% 1D, input observations, unequal occurrences
mu_a=0;
v_a=1;
obs_a=normrnd(mu_a,sqrt(v_a),[1000 1]);

mu_b=2;
v_b=1;
obs_b=normrnd(mu_b,sqrt(v_b),[2000 1]);

results=classify(obs_a,obs_b,'type','obs')

%% 2D, simple

mu_a=[4; 5];
v_a=[2 1.1; 1.1 1];

mu_b=[4; 5];
v_b=1.4*v_a;

results=classify([mu_a,v_a],[mu_b,v_b]);

%% 2D, input observations

mu_a=[2 4];
v_a=[1 1.5; 1.5 3];
obs_a = mvnrnd(mu_a,v_a,100);

mu_b=[5 0];
v_b=[3 0; 0 1];
obs_b = mvnrnd(mu_b,v_b,100);

results1=classify(obs_a,obs_b,'type','obs');

axis image
xlim([-10 10])
ylim([-10 10])

% modify boundary
custom_bd_coeffs=results1.bd_coeffs_obs_opt;
custom_bd_coeffs.a2=custom_bd_coeffs.a2+.2;
custom_bd_coeffs.a1=custom_bd_coeffs.a1-5;
custom_bd_coeffs.a0=custom_bd_coeffs.a0+10;

results2=classify(obs_a,obs_b,'type','obs','custom_bd_coeffs',custom_bd_coeffs);

axis image
xlim([-10 10])
ylim([-10 10])

%% 2D, camouflage

load('C:\Users\abhra\Google Drive\Geisler Lab\camouflage-detection\global_data\seed_energy_pn_gradbynorm.mat')
load('C:\Users\abhra\Google Drive\Geisler Lab\camouflage-detection\global_data\edge_spectra_LLR.mat')
load('C:\Users\abhra\Google Drive\Geisler Lab\camouflage-detection\global_data\seed_energy_pink_noise_whitened.mat')

n_obs=size(LLR_spectrum_pink_b,1);
obs_a = [seed_energy(1:n_obs,2),LLR_spectrum_pink_b];
obs_b = [seed_energy(1:n_obs,3),LLR_spectrum_pink_t];

results_pink=classify(obs_a,obs_b,'type','obs');

axis normal
xlim([-1 1])
ylim([-.5 3])

obs_a = [seed_energy_w(1:n_obs,2),LLR_spectrum_white_b];
obs_b = [seed_energy_w(1:n_obs,3),LLR_spectrum_white_t];

results_white=classify(obs_a,obs_b,'type','obs');

axis normal
xlim([-1 1])
ylim([-.5 3])

%% 3D, simple

mu_a=[0;0;0];
v_a=eye(3);

mu_b=[2;1;1];
v_b=2*eye(3);

results=classify([mu_a,v_a],[mu_b,v_b]);

%% 3D, complex

mu_a=[0;0;0];
v_a=[1 .5 0; .5 2 1; 0 1 4];

mu_b=[2; 1; 2];
v_b=[2 1 -1.5; 1 3 -2; -1.5 -2 4];

results=classify([mu_a,v_a],[mu_b,v_b]);

