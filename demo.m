% Examples for classify library.
% Author:
%   Abhranil Das <abhranil.das@utexas.edu>
%	Center for Perceptual Systems, University of Texas at Austin
% If you use this code, please cite:
%   A new method to compute classification error
%   https://jov.arvojournals.org/article.aspx?articleid=2750251

%% 1D, integrate 
mu=0;
v=1;

% integrate in a quadratic region (x-1)^2<1, i.e. f(x)= -1*x^2 + 1*x + 1 >0
reg_quad.a2=-1;
reg_quad.a1=1;
reg_quad.a0=1;
integrate_normal(mu,v,reg_quad);

% integrate in a region defined by a non-quadratic f(x)>0
reg_func=@(x) cos(x.^2);
integrate_normal(mu,v,reg_func,'reg_type','cheb','cheb_reg_span',10,'func_crossings',200);

%% 1D, classify
mu_1=0;
v_1=1;

mu_2=2.5;
v_2=1.5;

results=classify_normals([mu_1,v_1],[mu_2,v_2])

% with unequal priors
results=classify_normals([mu_1,v_1],[mu_2,v_2],'prior_1',.7)

% with outcome values
results=classify_normals([mu_1,v_1],[mu_2,v_2],'vals',[3 0; 0 1])

%% 1D, classify samples (priors are taken to be prop. to sample sizes)
mu_1=0;
v_1=1;
samp_1=normrnd(mu_1,sqrt(v_1),[700 1]);

mu_2=1.5;
v_2=1.5;
samp_2=normrnd(mu_2,sqrt(v_2),[300 1]);

results=classify_normals(samp_1,samp_2,'type','samp','opt_samp',false)

%% 2D, one inside the other
mu_1=[4; 5];
v_1=[2 1.1; 1.1 1];

mu_2=mu_1;
v_2=1.4*v_1;

results=classify_normals([mu_1,v_1],[mu_2,v_2])

%% 2D, classify samples
n_samp=1e4;

mu_1=[2 4];
v_1=[1 1.5; 1.5 3];
samp_1 = mvnrnd(mu_1,v_1,n_samp);

mu_2=[5 0];
v_2=[3 0; 0 1];
samp_2 = mvnrnd(mu_2,v_2,n_samp);

results=classify_normals(samp_1,samp_2,'type','samp')

axis image
xlim([-10 10])
ylim([-10 10])

% modify boundary
custom_reg_quad=results.samp_opt_reg_quad;
custom_reg_quad.a2=custom_reg_quad.a2+.1;
custom_reg_quad.a1=custom_reg_quad.a1-2.5;
custom_reg_quad.a0=custom_reg_quad.a0+5;

results=classify_normals(samp_1,samp_2,'type','samp','reg',custom_reg_quad)

axis image
xlim([-10 10])
ylim([-10 10])

%% 2D, classify non-normal samples
n_samp=1e3;
samp_1=exp(mvnrnd([0 0],eye(2),n_samp));
samp_2=exp(mvnrnd([1 1],eye(2),n_samp));

results=classify_normals(samp_1,samp_2,'type','samp')
xlim([-10 20])
ylim([-5 15])

%% Invert and combine regions of integration
mu=[0;0];
v=eye(2);

circle_left.a2=-eye(2);
circle_left.a1=[-2;0];
circle_left.a0=4;

circle_right=circle_left;
circle_right.a1=-circle_left.a1;

circle_left_rayscan=@(n,orig)ray_scan(circle_left,n,orig);
circle_right_rayscan=@(n,orig)ray_scan(circle_right,n,orig);

circle_union_rayscan=@(n,orig) combine_regs({circle_left_rayscan,circle_right_rayscan},'or',n,orig);
circle_intersection_rayscan=@(n,orig) combine_regs({circle_left_rayscan,circle_right_rayscan},'and',n,orig);

integrate_normal(mu,v,circle_union_rayscan,'reg_type','ray_scan');
axis image
xlim([-4 4])
ylim([-4 4])

figure
integrate_normal(mu,v,circle_intersection_rayscan,'reg_type','ray_scan');
axis image
xlim([-4 4])
ylim([-4 4])

circle_left_invert_rayscan=@(n,orig) invert_reg(circle_left_rayscan,n,orig);

crescent_rayscan=@(n,orig) combine_regs({circle_left_invert_rayscan,circle_right_rayscan},'and',n,orig);

figure
integrate_normal(mu,v,crescent_rayscan,'reg_type','ray_scan');
axis image
xlim([-4 4])
ylim([-4 4])

%% 2D, integrate in a region defined by f(x,y)>0
mu=[0;-1];
v=[1 .5; .5 2]/6;
reg_func=@(x,y) exp(-x.^2).*abs(cos(10*x))-y;
integrate_normal(mu,v,reg_func,'reg_type','cheb','cheb_reg_span',4,'RelTol',1e-1);
axis image
xlim([-2 2])
ylim([-2 2])

% combine two such regions
reg2_func=@(x,y) y+2*exp(-x.^2);

reg_rayscan=@(n,orig)ray_scan(reg_func,n,orig,'reg_type','cheb','cheb_reg_span',3);
reg2_rayscan=@(n,orig)ray_scan(reg2_func,n,orig,'reg_type','cheb','cheb_reg_span',3);
reg_rayscan_combined=@(n,orig) combine_regs({reg_rayscan,reg2_rayscan},'and',n,orig);

figure
integrate_normal(mu,v,reg_rayscan_combined,'reg_type','ray_scan');
axis image
xlim([-2 2])
ylim([-2.5 1.5])

%% 3D, simple
mu_1=[0;0;0];
v_1=eye(3);

mu_2=[2;1;1];
v_2=2*eye(3);

results=classify_normals([mu_1,v_1],[mu_2,v_2])

% plot the boundary points used for integration
hold on
plot3(results.norm_bd_pts(1,:),results.norm_bd_pts(2,:),results.norm_bd_pts(3,:),'.','markersize',4)

%% High-accuracy estimation of tiny errors (large d')
format long
dprime_true=75

mu_1=[0;0;0];
v_1=eye(3);

mu_2=dprime_true*[1;0;0];
v_2=(1+1e-12)*eye(3);

results=classify_normals([mu_1,v_1],[mu_2,v_2],'AbsTol',0,'RelTol',1e-3)
dprime_computed=results.norm_dprime
format

%% Classify data from actual detection experiment

dataArray = textscan(fopen('absent.txt','r'), '%*q%f%f%f%[^\n\r]', 'Delimiter', ',', 'HeaderLines' ,1);
absent = [dataArray{1:end-1}];

dataArray = textscan(fopen('present.txt','r'), '%*q%f%f%f%[^\n\r]', 'Delimiter', ',', 'HeaderLines' ,1);
present = [dataArray{1:end-1}];

results=classify_normals(absent,present,'type','samp')
axis normal
xlim([-1 2]); ylim([0 .5]); zlim([-400 400]);
xlabel('edge'); ylabel('luminance'); zlabel('pattern');

%% 3D, integrate in a region defined by f(x,y,z)>0
mu=[0;0;0];
v=eye(3);
reg_func=@(x,y,z) 3-(4-(x.^2+y.^2).^0.5).^2-z.^2;
[~,~,bd_pts]=integrate_normal(mu,v,reg_func,'reg_type','cheb','cheb_reg_span',20);
axis image

%% 4D classification
mu_1=[0;0;0;0];
v_1=eye(4);

mu_2=[1;1;1;1];
v_2=eye(4);

n_samp=1e3;
results=classify_normals(mvnrnd(mu_1,v_1,n_samp),...
    mvnrnd(mu_2,v_2,n_samp),'type','samp','prior_1',.7,'vals',[2 0; 0 1])

%% Integrate non-quadratic function f of a normal,
% equivalent to integrating normal in the non-quadratic region f>0, using chebfun

% first install chebfun
% unzip('https://github.com/chebfun/chebfun/archive/master.zip')
% movefile('chebfun-master', 'chebfun'), addpath(fullfile(cd,'chebfun')), savepath

mu_1=[2;3];
v_1=[1 -.5; -.5 2];
e=12;

e_cheb=@(x,y) e-x.*y.^2/2; % define f (as a cheb function)

% syms x y
% diff(e_cheb(x,y))
integrate_normal(mu_1,v_1,e_cheb,'reg_type','cheb');

% classify normals and samples wrt this region
mu_2=[3;4];
v_2=[2 -1; -1 4];
n_samp=1e3;
samp_1=mvnrnd(mu_1,v_1,n_samp);
samp_2=mvnrnd(mu_2,v_2,n_samp);
results=classify_normals(samp_1,samp_2,'prior_1',.8,'type','samp','reg',e_cheb,'reg_type','cheb');

% integrate multi-valued function of a normal
% integrate [f,g]=[e_cheb,p_cheb] where both are +ve
p=7;
p_quad.a2=(eye(2)-1)/2;
p_quad.a1=[0;0];
p_quad.a0=p;

e_rayscan=@(n,orig)ray_scan(e_cheb,n,orig,'reg_type','cheb','cheb_reg_span',15);
p_rayscan=@(n,orig)ray_scan(p_quad,n,orig);
ep_rayscan=@(n,orig) combine_regs({e_rayscan,p_rayscan},'and',n,orig);

figure; hold on
plot_boundary(e_cheb,2,'reg_type','cheb','cheb_reg_span',15,'plot_color','b')
plot_boundary(p_quad,2,'reg_type','quad','plot_color','r')
xlim([-5 15])
ylim([-15 15])

% classify normals and samples using this region
%results=classify_normals([mu_1,v_1],[mu_2,v_2],'prior_1',.8,'reg',{ep_rayscan,@(n,orig) invert_reg(ep_rayscan,n,orig)},'reg_type','ray_scan','n_bd_pts',1e3);
results=classify_normals(samp_1,samp_2,'prior_1',.8,'type','samp','reg',{ep_rayscan,@(n,orig) invert_reg(ep_rayscan,n,orig)},'reg_type','ray_scan','n_samp_bd_pts',1e3);
xlim([-5 20])
ylim([-15 15])

%% Multi-class classification, 1D
dists=struct;
dists(1).mu=-1; dists(1).v=1;
dists(2).mu=1; dists(2).v=2;
dists(3).mu=4; dists(3).v=.5;

% use priors
results=classify_normals_multi(dists,'priors',[.4 .5 .1])

% use outcome values
vals=[4 0 0; 0 1 -2; 0 0 1];
results=classify_normals_multi(dists,'priors',[.4 .5 .1],'vals',vals)
%% Multi-class, 2D

% define means and vcovs of the normals
mus=[[1;0],[0;1],[-1;0],[0;-1]];
vs=cat(3,2*eye(2),eye(2),eye(2),eye(2));

% define struct of all normals
dists=struct;
for i=1:4
    dists(i).mu=mus(:,i); dists(i).v=vs(:,:,i);
end

% plot the multi-class boundary of normal 3
plot_boundary(@(n,orig) opt_reg_multi(n,mus,vs,'idx',3,'orig',orig),2,'orig',mus(:,3),'reg_type','ray_scan')
title 'boundary of normal 3'

% classify
results=classify_normals_multi(dists)
axis image

% now input samples from these normals
dists2=struct;
for i=1:4
    dists2(i).sample=mvnrnd(dists(i).mu,dists(i).v,1e2);
end
results_samp=classify_normals_multi(dists2,'type','samp')
axis image

%% Multi-class, 3D
dists=struct;
dists(1).mu=[1;0;0]; dists(1).v=2*eye(3);
dists(2).mu=[0;1;0]; dists(2).v=eye(3);
dists(3).mu=[-1;0;0]; dists(3).v=eye(3);
dists(4).mu=[0;-1;0]; dists(4).v=eye(3);

results=classify_normals_multi(dists)
