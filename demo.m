% Examples for classify library.
% Credits:
%   Abhranil Das <abhranil.das@utexas.edu>
%	R Calen Walshe
%	Wilson S Geisler
%	Center for Perceptual Systems, University of Texas at Austin
% If you use this code, please cite:
%   A new method to compute classification error
%   https://jov.arvojournals.org/article.aspx?articleid=2750251

%% 1D, simple
mu_a=0;
v_a=1;

mu_b=2.5;
v_b=1.5;

results=classify_normals([mu_a,v_a],[mu_b,v_b]);

% with outcome utilities
results=classify_normals([mu_a,v_a],[mu_b,v_b],'val',[2 0; 0 1]);

%% 1D, unequal prior
mu_a=0;
v_a=1;

mu_b=0.5;
v_b=1.5;

results=classify_normals([mu_a,v_a],[mu_b,v_b],'p_a',.7);

%% 1D, input observations, unequal occurrences
mu_a=0;
v_a=1;
obs_a=normrnd(mu_a,sqrt(v_a),[700 1]);

mu_b=0.5;
v_b=1.5;
obs_b=normrnd(mu_b,sqrt(v_b),[300 1]);

results=classify_normals(obs_a,obs_b,'type','obs');

%% 1D, input observations, outcome values
mu_a=0;
v_a=1;
obs_a=normrnd(mu_a,sqrt(v_a),[1000 1]);

mu_b=2.5;
v_b=1.5;
obs_b=normrnd(mu_b,sqrt(v_b),[1000 1]);

results=classify_normals(obs_a,obs_b,'type','obs',...
    'val',[4 0; 0 1]);

%% 2D, simple
mu_a=[4; 5];
v_a=[2 1.1; 1.1 1];

mu_b=mu_a;
v_b=1.4*v_a;

results=classify_normals([mu_a,v_a],[mu_b,v_b]);

%% 2D, input observations
n_samp=1e2;

mu_a=[2 4];
v_a=[1 1.5; 1.5 3];
obs_a = mvnrnd(mu_a,v_a,n_samp);

mu_b=[5 0];
v_b=[3 0; 0 1];
obs_b = mvnrnd(mu_b,v_b,n_samp);

results1=classify_normals(obs_a,obs_b,'type','obs');

axis image
xlim([-10 10])
ylim([-10 10])

% modify boundary
custom_bd_coeffs=results1.bd_coeffs_obs;
custom_bd_coeffs.a2=custom_bd_coeffs.a2+.2;
custom_bd_coeffs.a1=custom_bd_coeffs.a1-5;
custom_bd_coeffs.a0=custom_bd_coeffs.a0+10;

results2=classify_normals(obs_a,obs_b,'type','obs',...
    'custom_bd_coeffs',custom_bd_coeffs);

axis image
xlim([-10 10])
ylim([-10 10])

%% 2D, input non-Gaussian observations
n_samp=1e3;
obs_a=unifrnd(0,3,[n_samp 2]);
obs_b=unifrnd(1,6,[n_samp 2]);

results=classify_normals(obs_a,obs_b,'type','obs');
xlim([-1 5]); ylim([-1 5]); axis image

%% 3D, simple
mu_a=[0;0;0];
v_a=eye(3);

mu_b=[2;1;1];
v_b=2*eye(3);

tic
results=classify_normals([mu_a,v_a],[mu_b,v_b]);
toc

% force grid method
bd_fn_a=@(n) opt_bd(n,[mu_a,v_a],[mu_b,v_b]);
bd_fn_b=@(n) opt_bd(n,[mu_b,v_b],[mu_a,v_a]);

tic
results_grid=classify_normals([mu_a,v_a],[mu_b,v_b],'custom_bd_fns',{bd_fn_a,bd_fn_b});
toc

%% 3D, far apart
mu_a=[0;0;0];
v_a=eye(3);

mu_b=[100;0;0];
v_b=1.5*eye(3);

results=classify_normals([mu_a,v_a],[mu_b,v_b]);

%% 3D, complex
mu_a=[0;0;0];
v_a=[1 .5 0; .5 2 1; 0 1 4];

mu_b=[2; 1; 2];
v_b=[2 1 -1.5; 1 3 -2; -1.5 -2 4];

results=classify_normals([mu_a,v_a],[mu_b,v_b]);

%% 3D, from actual detection experiment

dataArray = textscan(fopen('absent.txt','r'), '%*q%f%f%f%[^\n\r]', 'Delimiter', ',', 'HeaderLines' ,1);
absent = [dataArray{1:end-1}];

dataArray = textscan(fopen('present.txt','r'), '%*q%f%f%f%[^\n\r]', 'Delimiter', ',', 'HeaderLines' ,1);
present = [dataArray{1:end-1}];

results=classify_normals(absent,present,'type','obs');
axis normal
xlim([0 1.5]); ylim([0 .5]); zlim([-200 200]);
xlabel('edge'); ylabel('luminance'); zlabel('pattern');

%% 4D, simple
mu_a=[0;0;0;0];
v_a=eye(4);

mu_b=[1;1;1;1];
v_b=2*eye(4);

results=classify_normals([mu_a,v_a],[mu_b,v_b]);

%% 4D, input observations
mu_a=[0;0;0;0];
v_a=eye(4);
obs_a = mvnrnd(mu_a,v_a,1000);

mu_b=[1;1;1;1];
v_b=eye(4);
obs_b = mvnrnd(mu_b,v_b,1000);

results=classify_normals(obs_a,obs_b,'type','obs');

