% Examples for classify library.
% Author:
%   Abhranil Das <abhranil.das@utexas.edu>
%	Center for Perceptual Systems, University of Texas at Austin
% If you use this code, please cite:
%   A new method to compute classification error
%   https://jov.arvojournals.org/article.aspx?articleid=2750251

%% 1D, integrate 
mu=-2;
v=3;

% integrate in a quadratic domain x^2-x-1<0, i.e. f(x)= -1*x^2 + 1*x + 1 >0
dom_quad.q2=-1;
dom_quad.q1=1;
dom_quad.q0=1;
integrate_normal(mu,v,dom_quad);

% most plots can be zoomed and panned.

% integrate in a domain defined by a non-quadratic f(x)>0
fun=@(x) cos(x.^2);
figure;
integrate_normal(mu,v,fun,'dom_type','fun','fun_span',5,'fun_resol',500);

% plot the pdf and cdf of f(x)
x=linspace(-2,2,100);
f=norm_fun_pdf(x,mu,v,fun,'fun_span',5,'fun_resol',500,'dx',1e-2);
figure; plot(x,f);
F=norm_fun_cdf(x,mu,v,fun,'fun_span',5,'fun_resol',500);
figure; plot(x,F);

%% 1D, classify
mu_1=0;
v_1=1;

mu_2=2.5;
v_2=1.5;

results=classify_normals([mu_1,v_1],[mu_2,v_2])

% with unequal priors
results=classify_normals([mu_1,v_1],[mu_2,v_2],'prior_1',.7)

% with outcome values
results=classify_normals([mu_1,v_1],[mu_2,v_2],'prior_1',.7,'vals',[3 0; 0 1])

%% 1D, classify using samples (priors are assumed prop. to sample sizes)
mu_1=0;
v_1=1;
samp_1=normrnd(mu_1,sqrt(v_1),[700 1]);

mu_2=1.5;
v_2=1.5;
samp_2=normrnd(mu_2,sqrt(v_2),[300 1]);

results=classify_normals(samp_1,samp_2,'type','samp')

%% 2D, integrate
mu=[-1; -1];
v=[1 0.5; 0.5 2];

% integrate in a quadratic domain (x+y)^2 > x+1, i.e.
% q(x,y)= [x y]*[1 1; 1 1]*[x;y] + [-1 0]*[x;y] -1 > 0
dom_quad.q2=[1 1; 1 1];
dom_quad.q1=[-1;0];
dom_quad.q0=-1;

% compare two integration algorithms
figure; integrate_normal(mu,v,dom_quad); % ray method
xlim([-6 10]); ylim([-10 2])
figure; integrate_normal(mu,v,dom_quad,'method','gx2'); % gx2 method
xlim([-6 10]); ylim([-10 2])

%% 2D, classify two normals, one inside the other
mu_1=[4; 5];
v_1=[2 1; 1 1];

mu_2=mu_1;
v_2=3*[2 -1; -1 1];

% compare two integration algorithms
results_ray=classify_normals([mu_1,v_1],[mu_2,v_2]) % ray method
xlim([-2 8]); ylim([-0 10]);
results_gx2=classify_normals([mu_1,v_1],[mu_2,v_2],'method','gx2') % gx2 method
xlim([-2 8]); ylim([-0 10]);

%% PAPER 2D, classify with custom boundaries and from samples

mu_1=[2;4];
v_1=[1 1.5; 1.5 3];

mu_2=[5;0];
v_2=[3 0; 0 1];

% ray method
results_ray=classify_normals([mu_1,v_1],[mu_2,v_2]); 
axis image; xlim([-10 10]); ylim([-10 10])
results_ray.norm_err

% compare with generalized chi square method
results_gx2=classify_normals([mu_1,v_1],[mu_2,v_2],'method','gx2'); 
axis image; xlim([-10 10]); ylim([-10 10])
results_gx2.norm_err

% now supply a custom linear boundary
linear_bd.q2=zeros(2);
linear_bd.q1=[-.7;1];
linear_bd.q0=0;

hold on;
plot_boundary(linear_bd,2,'dom_type','quad','mu',[-4;4],'plot_type','line','line_color',[1 0 1])
axis image; xlim([-10 10]); ylim([-10 10])
set(gca,'fontsize',13); box off

results_linear=classify_normals([mu_1,v_1],[mu_2,v_2],'dom',linear_bd,'plotmode',false)

% now classify using samples
n_samp=1e3;
samp_1=mvnrnd(mu_1,v_1,n_samp);
samp_2=mvnrnd(mu_2,v_2,n_samp);

results_samp=classify_normals(samp_1,samp_2,'type','samp')
axis image; xlim([-10 10]); ylim([-10 10])

% modify the sample-optimized boundary
custom_bd=results_samp.samp_opt_bd;
custom_bd.q2=custom_bd.q2+.1;
custom_bd.q1=custom_bd.q1-2.5;
custom_bd.q0=custom_bd.q0+5;

results_samp_custom=classify_normals(samp_1,samp_2,'type','samp','dom',custom_bd)
axis image; xlim([-10 10]); ylim([-10 10])

%% PAPER 2D, classify non-normal samples
n_samp=1e3;
samp_1=exp(mvnrnd([0 0],eye(2),n_samp));
samp_2=-exp(mvnrnd([1 1],eye(2),n_samp));

results=classify_normals(samp_1,samp_2,'type','samp')
axis image; xlim([-20 15]); ylim([-20 15])
set(gca,'fontsize',13); box off

%% Inversion/union/intersection of integration/classification domains
mu=[0;0];
v=[.5 0; 0 1];

circle_left=@(x,y) -(x+1).^2-y.^2+5;
circle_right=@(x,y) -(x-1).^2-y.^2+5;

circle_union=@(x,y) max(circle_left(x,y), circle_right(x,y));
integrate_normal(mu,v,circle_union,'dom_type','fun','fun_span',5);
axis image; xlim([-4 4]); ylim([-4 4])

circle_intersection=@(x,y) min(circle_left(x,y), circle_right(x,y));
figure
integrate_normal(mu,v,circle_intersection,'dom_type','fun','fun_span',5);
axis image; xlim([-4 4]); ylim([-4 4])

crescent=@(x,y) max(circle_left(x,y), -circle_right(x,y));
figure
integrate_normal(mu,v,crescent,'dom_type','fun','fun_span',5);
axis image; xlim([-4 4]); ylim([-4 4])

% classify normals using this domain
mu_2=[2.2;0];
v_2=[.5 0; 0 .25];
classify_normals([mu,v],[mu_2,v_2],'dom',crescent,'dom_type','fun','fun_span',5);
axis image; xlim([-2 4]); ylim([-3 3])
set(gca,'fontsize',13); box off

%% 3D, classify
mu_1=[0;0;0];
v_1=eye(3);

mu_2=[2;1;1];
v_2=2*eye(3);

results=classify_normals([mu_1,v_1],[mu_2,v_2])

%% High-accuracy estimation of tiny errors (large d')
format long
dprime_true=75

mu_1=[0;0;0];
v=eye(3);

mu_2=dprime_true*[1;0;0];

results=classify_normals([mu_1,v],[mu_2,v],'method','gx2','AbsTol',0,'RelTol',1e-2)
axis image; xlim([-5 80]); ylim([-2 2]); zlim([-2 2]); view(-53,25)
dprime_computed=results.norm_dprime
format

%% PAPER accuracy vs separation between the normals

mu_1=[0;0;0];

% both normals have the same covariance, so we can check against
% the true d' (Mahalanobis distance).
v=[1 .5 .7;
  .5  2  1 ;
  .7  1  3];

steps=linspace(1,100,20);
d_true=nan(size(steps));
d_gx2=nan(size(steps));
d_ray=nan(size(steps));

parfor i=1:length(steps)
    i
    mu_2=steps(i)*[1;1;1];
    
    results_gx2=classify_normals([mu_1,v],[mu_2,v],'method','gx2','AbsTol',0,'RelTol',0,'plotmode',false);
    d_true(i)=results_gx2.norm_maha_dprime;
    d_gx2(i)=results_gx2.norm_dprime;
    
    results_ray=classify_normals([mu_1,v],[mu_2,v],'method','ray','AbsTol',0,'RelTol',0,'plotmode',false);
    d_ray(i)=results_ray.norm_dprime;
end

rel_err_gx2=abs(d_gx2-d_true)./d_true;
rel_err_ray=abs(d_ray-d_true)./d_true;

figure; hold on
plot(d_true,rel_err_gx2,'-o')
plot(d_true,rel_err_ray,'-o')
xline(-2*norminv(realmin)) % largest computable d', corr. to the smallest possible error representable in double-precision 
yline(eps) % machine epsilon for double precision
xlim([1 80]); ylim([0 2.2e-15])
xlabel 'true d'''
ylabel('$\epsilon$','interpreter','latex')
legend({'$\tilde{\chi}^2$','ray'},'interpreter','latex');
set(gca,'fontsize',13); box off; legend boxoff

%% PAPER Integrate in a toroidal domain defined by implicit function f(x)>0

mu=[0;0;0];
v=[1 0 0;
   0 8 4;
   0 4 8];

% plot the error ellipsoid of the normal
plot_normal(mu,v);

fun_torus=@(x1,x2,x3) 1.5-(5-(x1.^2+x2.^2).^0.5).^2-x3.^2;
figure;
integrate_normal(mu,v,fun_torus,'dom_type','fun','fun_span',3,'fun_resol',10,'RelTol',1e-1);
axis image; xlim([-7 7]); ylim([-7 7]); zlim([-7 7]);
set(gca,'xtick',linspace(-7,7,5)); set(gca,'ytick',linspace(-7,7,5)); set(gca,'ztick',linspace(-7,7,5))
set(gca,'fontsize',13)

% plot the pdf and cdf of f(x)
x=[linspace(-80,-20,30) linspace(-20,10,20)];
pdf=norm_fun_pdf(x,mu,v,fun_torus,'fun_span',3,'fun_resol',10,'RelTol',1e-1,'dx',8);
colors=colororder;
figure; area(x,pdf,'facecolor',colors(1,:),'facealpha',0.4,'edgecolor',colors(1,:))
set(gca,'ytick',[]); set(gca,'fontsize',13); box off
cdf=norm_fun_cdf(x,mu,v,fun_torus,'fun_span',3,'fun_resol',10,'RelTol',1e-1);
figure; plot(x,cdf)

%% PAPER 4D, integrate
mu=[1;1;1;1];
v=diag([1 2 3 4]);

% integrate outside the sphere of radius 5, i.e.
% x1^2+x2^2+x3^2+x4^2-5^2>0, i.e. x'*eye(4)*x + zeros(4,1)'*x -25 >0
dom_quad.q2=eye(4);
dom_quad.q1=zeros(4,1);
dom_quad.q0=-25;

p=integrate_normal(mu,v,dom_quad)
xlim([-30 40]); ylim([0 .06]);
set(gca,'ytick',[])
set(gca,'fontsize',13); box off

%% PAPER 4D, classify, with priors and outcome values
mu_1=[0;0;0;0];
v_1=eye(4);

mu_2=[1;1;1;1];
v_2=eye(4);

n_samp=1e3;

% here the Bayes decision variable is normally distributed
results=classify_normals([mu_1,v_1],[mu_2,v_2],'prior_1',.7,'vals',[2 0; 0 1])

mu_1=[0;0;0;0];
v_1=[1 0 0 0;
     0 2 -1 0;
     0 -1 3 2;
     0 0 2 4];

mu_2=1.5*[1;1;1;1];
v_2=[2 0 0 0;
     0 2 0 0;
     0 0 2 0;
     0 0 0 1];

% correct classification of class 1 is valued 4x than class 2
% here the Bayes decision variable is distributed as a generalized chi-square
results=classify_normals([mu_1,v_1],[mu_2,v_2],'prior_1',.7,'vals',[4 0; 0 1])
results=classify_normals([mu_1,v_1],[mu_2,v_2],'prior_1',.7,'vals',[4 0; 0 1],'method','ray','mc_samples',5e2,'plotmode',false)
ylim([0 .12])

% now classify using samples
n_samp=1e3;
results=classify_normals(mvnrnd(mu_1',v_1,n_samp),mvnrnd(mu_2',v_2,n_samp),'type','samp','prior_1',.7,'vals',[4 0; 0 1])
set(gca,'fontsize',13); box off

%% PAPER 4D Integrate in an octahedral domain, using Monte Carlo
mu=zeros(4,1);
v=eye(4);

fun=@(x1,x2,x3,x4) 1-abs(x1)-abs(x2)-abs(x3)-abs(x4);
mc_samples=round(10.^linspace(1,4,20)); % # of Monte Carlo samples
n_repeat=5;
plist=nan(length(mc_samples),n_repeat);

for k=1:n_repeat
    for i=1:length(mc_samples)
        [k i]
        p=integrate_normal(mu,v,fun,'dom_type','fun','fun_span',3,'fun_resol',10,'mc_samples',mc_samples(i))
        plist(i,k)=p;
    end
end

plot(mc_samples,plist,'-k','marker','.','markersize',10)
xlabel 'Monte Carlo sample size'
ylabel 'p'
set(gca,'xscale','log')
set(gca,'fontsize',13); box off

%% PAPER Integrate vector-valued function of a normal

% (x,y) is a normal vector with these parameters:
mu=[6;-19];
v=[1 -.7; -.7 2];

% functions of the normal vector
f1=@(x,y) cos(x);
f2=@(x,y) cos(y);

% integrate f1 above 0
figure;
integrate_normal(mu,v,f1,'dom_type','fun');
axis image; xlim([-10 20]); ylim([-50 10]);

% integrate f2 above 0
figure;
integrate_normal(mu,v,f2,'dom_type','fun');
axis image; xlim([-10 20]); ylim([-50 10]);

% prob. that cos(x) and cos(y) are both >0
% i.e. integrate vector function f=[f1,f2] in the domain f1>0 and f2>0, i.e. min(f1,f2)>0.
f_domain1=@(x,y) min(f1(x,y),f2(x,y));

figure;
integrate_normal(mu,v,f_domain1,'dom_type','fun','fun_span',10);
axis image; xlim([-10 20]); ylim([-50 10]);
set(gca,'fontsize',13); box off

% prob. that cos(x) + cos(y) > 0.5
% i.e. integrate the vector function in the implicit domain f1+f2-0.5 > 0
f_domain2=@(x,y) f1(x,y)+f2(x,y)-.5;

figure;
integrate_normal(mu,v,f_domain2,'dom_type','fun','fun_span',10);
axis image; xlim([-10 20]); ylim([-50 10]);
set(gca,'fontsize',13); box off

%% Classifying 3 normals, 1D
normals=struct;
normals(1).mu=-1; normals(1).v=1;
normals(2).mu=1; normals(2).v=2;
normals(3).mu=4; normals(3).v=.5;

% priors and outcome values
priors=[.4 .5 .1];
vals=[4 0 0; 0 1 -2; 0 0 1];

results=classify_normals_multi(normals,'priors',priors,'vals',vals)

%% Classifying 4 normals, 2D

% define struct of all normals
normals=struct;
normals(1).mu=[1;0]; normals(1).v=2*eye(2);
normals(2).mu=[0;1]; normals(2).v=eye(2);
normals(3).mu=[-1;0]; normals(3).v=eye(2);
normals(4).mu=[0;-1]; normals(4).v=eye(2);

% plot the multi-class boundary of normal 3
plot_boundary(@(n,mu,v) opt_class_multi(n,normals,3,'mu',mu,'v',v),2,'mu',normals(3).mu,'dom_type','ray_scan')
title 'boundary of normal 3'

% classify
results=classify_normals_multi(normals)
axis image

q12=opt_class_quad([normals(1).mu normals(1).v],[normals(2).mu normals(2).v]);
q13=opt_class_quad([normals(1).mu normals(1).v],[normals(3).mu normals(3).v]);
f12=quad2fun(q12,1);
f13=quad2fun(q13,1);
f=@(x,y) min(f12(x,y),f13(x,y));
plot_boundary(f,2,'dom_type','fun','plot_type','line')

% now classify using samples from these normals
samples=struct;
for i=1:4
    samples(i).sample=mvnrnd(normals(i).mu,normals(i).v,1e4);
end
results_samp=classify_normals_multi(samples,'type','samp')
axis image

% now classify using samples from t distributions similar to these normals
samples=struct;
for i=1:4
    samples(i).sample=mvtrnd(normals(i).v,3,1e4)+normals(i).mu';
end
results_samp=classify_normals_multi(samples,'type','samp')
axis image; xlim([-6 6]); ylim([-6 6]);

%% PAPER Classifying 7 normals, 2D

% define struct of all normals
normals=struct;
normals(1).mu=[2;0]; normals(1).v=[2 1; 1 2];
normals(2).mu=[0;1]; normals(2).v=[.5 0; 0 1];
normals(3).mu=[-1;0]; normals(3).v=[1 .3; .3 1];
normals(4).mu=[0;-1]; normals(4).v=[.5 -.5; -.5 1];
normals(5).mu=[-2;2]; normals(5).v=.5*[1 1; 1 5];
normals(6).mu=[2;-3]; normals(6).v=.3*eye(2);
normals(7).mu=[-2;-2.5]; normals(7).v=[1 0; 0 .1];

% classify
results=classify_normals_multi(normals)
axis image; xlim([-10 10]); ylim([-10 10]);
set(gca,'fontsize',13); box off

%% Classifying 4 normals, 3D
normals=struct;
normals(1).mu=[1;0;0]; normals(1).v=2*eye(3);
normals(2).mu=[0;1;0]; normals(2).v=eye(3);
normals(3).mu=[-1;0;0]; normals(3).v=eye(3);
normals(4).mu=[0;-1;0]; normals(4).v=eye(3);

results=classify_normals_multi(normals)

%% PAPER Classifying 4 normals, 4D
% normals=struct;
% normals(1).mu=[0;0;-1;-1]; normals(1).v=diag([1 2 3 4]);
% normals(2).mu=[1;1;0;0]; normals(2).v=eye(4)/4;
% normals(3).mu=[2;2;1;1]; normals(3).v=eye(4);
% normals(4).mu=[3;3;3;3]; normals(4).v=eye(4);
% 
% priors=[.3 .3 .35 .05];
% results=classify_normals_multi(normals,'priors',priors,'mc_samples',1e3,'plotmode',[1;1;1;1])
% 
% % now classify using samples from these normals
% samples=struct;
% for i=1:4
%     samples(i).sample=mvnrnd(normals(i).mu,normals(i).v,2e2);
% end
% results_samp=classify_normals_multi(samples,'type','samp','priors',priors,'mc_samples',1e3,'plotmode',[1;1;1;1])
% set(gca,'fontsize',13); box off

%% PAPER Classifying 4 4D t distribution samples
params=struct;
params(1).mu=[0;-4;-2;-2]; params(1).v=diag([1 2 3 4]);
params(2).mu=[1;1;0;0]; params(2).v=eye(4)/4;
params(3).mu=[3;4;5;7]; params(3).v=eye(4);
params(4).mu=[5;5;7;4]; params(4).v=eye(4);

priors=[.3 .5 .15 .05];

n_samp=1e3;
samples=struct;
for i=1:4
    samples(i).sample=mvtrnd(params(i).v,3,n_samp)+params(i).mu';
end
results_samp=classify_normals_multi(samples,'type','samp','priors',priors,'mc_samples',1e3,'plotmode',[1;1;1;1])
set(gca,'fontsize',13); box off

%% PAPER testing the normal approximation for 2 2d classes
N_samp=1e4; n_samp=1e2;

mu_1=[5 5]; v_1=[1 .5; .5 1]; samp_1=mvnrnd(mu_1,v_1,N_samp).^3; % not normal
mu_2=[-200 -200]; v_2=[1 1; 1 3]*2e3; samp_2=mvnrnd(mu_2,v_2,N_samp);

results=classify_normals(samp_1,samp_2,'type','samp','samp_opt',false);
xlim([-300 600]); ylim([-400 600]); hold on
bd=results.norm_bd;

q0s=linspace(-100,35,100)';
true_p11=nan(size(q0s)); true_p2_2=nan(size(q0s)); % p(1|1)
norm_p11=nan(size(q0s)); norm_p2_2=nan(size(q0s));

for i=1:length(q0s)
    i
    bd_shift=bd; bd_shift.q0=q0s(i);
    if ~rem(i,10)
        plot_boundary(bd_shift,2,'plot_type','line','line_color',.5*[1 1 1]);
    end
    results_shift=classify_normals(samp_1,samp_2,'type','samp','dom',bd_shift,'samp_opt',false,'plotmode',false);
    
    % true outcomes
    true_p11(i)=results_shift.samp_errmat(1,1)/N_samp;
    true_p2_2(i)=results_shift.samp_errmat(2,2)/N_samp;
    
    % normal outcomes
    norm_p11(i)=results_shift.norm_errmat(1,1)/.5;
    norm_p2_2(i)=results_shift.norm_errmat(2,2)/.5;
    
end
true_p1_1_sd=sqrt(true_p11.*(1-true_p11)/n_samp);
true_p2_2_sd=sqrt(true_p2_2.*(1-true_p2_2)/n_samp);

norm_p1_1_sd=sqrt(norm_p11.*(1-norm_p11)/n_samp);
norm_p2_2_sd=sqrt(norm_p2_2.*(1-norm_p2_2)/n_samp);

figure; hold on
colors=colororder;
xline(bd.q0,'k','optimal boundary','LabelVerticalAlignment','bottom') % optimal criterion

% p(1|1)
x=q0s; y=true_p11; dy=true_p1_1_sd;
fill([x;flipud(x)],[y-dy;flipud(y+dy)],colors(1,:),'facealpha',.3,'linestyle','none');

y=norm_p11; dy=norm_p1_1_sd;
fill([x;flipud(x)],[y-dy;flipud(y+dy)],'k','facecolor','none','edgecolor',colors(1,:));

% p(2|2)
y=true_p2_2; dy=true_p2_2_sd;
fill([x;flipud(x)],[y-dy;flipud(y+dy)],colors(2,:),'facealpha',.3,'linestyle','none');

y=norm_p2_2; dy=norm_p2_2_sd;
fill([x;flipud(x)],[y-dy;flipud(y+dy)],'k','facecolor','none','edgecolor',colors(2,:));

xlim([min(q0s) max(q0s)]); ylim([-.01 1.01]); xlabel('boundary offset q_0')
set(gca,'fontsize',13); box off

%% PAPER testing the normal approximation for 4 4d classes
n_class=4; N_samp=1e4; n_samp=1e2;

% t distribution parameters
params_t=struct;
params_t(1).mu=[0;-4;-2;-2]; params_t(1).v=diag([1 2 3 4]);
params_t(2).mu=[1;1;0;0]; params_t(2).v=eye(4)/4;
params_t(3).mu=[3;4;5;7]; params_t(3).v=eye(4);
params_t(4).mu=[5;5;7;4]; params_t(4).v=eye(4);

samples=struct; params_samp=struct;
for i=1:4
    % generate t samples
    sample=mvtrnd(params_t(i).v,3,N_samp)+params_t(i).mu';
    samples(i).sample=sample;
    % mean and covariance of samples
    params_samp(i).mu=mean(sample)';
    params_samp(i).v=cov(sample)';
end

len_fam=100;
prior_fam=10.^linspace(-6,12,len_fam)';
vscale_fam=10.^linspace(-1,3,len_fam)';

true_p22_mean=nan(len_fam,1);
true_pe_mean=nan(len_fam,1);
true_pe_sd=nan(len_fam,1);

norm_p22_mean=nan(len_fam,1);
norm_pe_mean=nan(len_fam,1);
norm_pe_sd=nan(len_fam,1);

parfor i=1:len_fam
    i
    % generate family of classification domains
    
    % uncomment when varying prior, comment when varying variance
    priors_shift=[1 prior_fam(i) 1 1]; priors_shift=priors_shift/sum(priors_shift);
    domains_shift=cell(n_class,1);
    for k=1:n_class
        domains_shift{k}=@(n,mu,v) opt_class_multi(n,params_samp,k,'mu',mu,'v',v,'priors',priors_shift);
    end
    
    % uncomment when varying variance, comment when varying prior
%     params_samp_shift=params_samp;
%     for k=1:n_class
%         params_samp_shift(k).v=vscale_fam(i)*params_samp(k).v;
%     end    
%     domains_shift=cell(n_class,1);
%     for k=1:n_class
%         domains_shift{k}=@(n,mu,v) opt_class_multi(n,params_samp_shift,k,'mu',mu,'v',v);
%     end
    
    results_shift=classify_normals_multi(samples,'type','samp','doms',domains_shift,'mc_samples',5e3,'plotmode',false);
    
    % true error rates
    true_p22_mean(i)=results_shift.samp_errmat(2,2)/N_samp;
    true_pe_mean(i)=results_shift.samp_err;
    samp_errmat_p=results_shift.samp_errmat/(n_class*N_samp);
    errs=sum(~eye(length(samples)).*samp_errmat_p,2);
    accs=diag(samp_errmat_p);
    true_pe_sd(i)=sqrt(sum(errs.*accs)/n_samp);
    
    % normal error rates
    norm_p22_mean(i)=results_shift.norm_errmat(2,2)/sum(results_shift.norm_errmat(2,:));
    norm_pe_mean(i)=results_shift.norm_err;
    errs=sum(~eye(length(samples)).*results_shift.norm_errmat,2);
    accs=diag(results_shift.norm_errmat);
    norm_pe_sd(i)=sqrt(sum(errs.*accs)/n_samp);
    
end
true_p22_sd=sqrt(true_p22_mean.*(1-true_p22_mean)/n_samp);
norm_p22_sd=sqrt(norm_p22_mean.*(1-norm_p22_mean)/n_samp);

figure; hold on;

% uncomment when varying prior, comment when varying variance
xline(1/3,'k','optimal boundary','LabelVerticalAlignment','bottom')
xlabel('assumed prior ratio p_1/(p_2+p_3+p_4)')
x=prior_fam/3;

% uncomment when varying variance, comment when varying prior
% xline(1,'k','optimal boundary','LabelVerticalAlignment','bottom')
% xlabel('assumed variance scale')
% x=vscale_fam;

% p(2|2)
y=true_p22_mean; dy=true_p22_sd;
fill([x;flipud(x)],[y-dy;flipud(y+dy)],[0.9290 0.6940 0.1250],'facealpha',.3,'linestyle','none');

y=norm_p22_mean; dy=norm_p22_sd;
fill([x;flipud(x)],[y-dy;flipud(y+dy)],'k','facecolor','none','edgecolor',[0.9290 0.6940 0.1250]);

% pe
y=true_pe_mean; dy=true_pe_sd;
fill([x;flipud(x)],[y-dy;flipud(y+dy)],[0.6 0.4 0.2],'facealpha',.3,'linestyle','none');

y=norm_pe_mean; dy=norm_pe_sd;
fill([x;flipud(x)],[y-dy;flipud(y+dy)],'k','facecolor','none','edgecolor',[0.6 0.4 0.2]);

xlim([min(x) max(x)]); set(gca,'xscale','log'); ylim([-.01 1.01])
set(gca,'fontsize',13); box off

%% PAPER Actual vision research data: detecting targets on natural scenes

absent=importdata('target_absent.txt',',',1);
present=importdata('target_present.txt',',',1);

results=classify_normals(absent.data,present.data,'type','samp')

axis normal
xlim([-50 170]); ylim([-320 200]); zlim([0 1000]); view(166,40);
xlabel('template'); ylabel('silhouette'); zlabel('edge');
set(gca,'fontsize',13); box off

%% PAPER Actual vision research data: detecting camouflage

load camouflage_edge_data
results_joint_2=classify_normals([edge_powers_2(:,1),edge_lpr_2(:,1)],[edge_powers_2(:,2),edge_lpr_2(:,2)],'type','samp');
xlim([0.2 1]); ylim([-900 400])
xlabel 'edge power'; ylabel 'edge spectrum'
set(gca,'fontsize',13); box off

results_dv_2=classify_normals(results_joint_2.samp_opt_dv{1},results_joint_2.samp_opt_dv{2},'type','samp','dom',@(x) x,'dom_type','fun','samp_opt',false);
xlim([-50 100]); ylim([0 .04]); set(gca,'ytick',[]);
set(gca,'fontsize',13); box off
xlabel('$q_s(${\boldmath$x$}$)$','interpreter','latex');

results_joint_4=classify_normals([edge_powers_4(:,1),edge_lpr_4(:,1)],[edge_powers_4(:,2),edge_lpr_4(:,2)],'type','samp','plotmode',false);
results_joint_8=classify_normals([edge_powers_8(:,1),edge_lpr_8(:,1)],[edge_powers_8(:,2),edge_lpr_8(:,2)],'type','samp','plotmode',false);

results_dv_joint=classify_normals([results_joint_2.samp_opt_dv{1},results_joint_4.samp_opt_dv{1},results_joint_8.samp_opt_dv{1}],[results_joint_2.samp_opt_dv{2},results_joint_4.samp_opt_dv{2},results_joint_8.samp_opt_dv{2}],'type','samp');
axis image; xlim([-100 125]); ylim([-50 75]); zlim([-50 50]); view(-20,24)
xlabel('$q_s(${\boldmath$x$}$)$ (2px)','interpreter','latex');
ylabel('$q_s(${\boldmath$x$}$)$ (4px)','interpreter','latex');
zlabel('$q_s(${\boldmath$x$}$)$ (8px)','interpreter','latex');
set(gca,'fontsize',13); box off

results_all=classify_normals([edge_powers_2(:,1),edge_lpr_2(:,1),...
                              edge_powers_4(:,1),edge_lpr_4(:,1),...
                              edge_powers_8(:,1),edge_lpr_8(:,1)],...
                             [edge_powers_2(:,2),edge_lpr_2(:,2),...
                              edge_powers_4(:,2),edge_lpr_4(:,2),...
                              edge_powers_8(:,2),edge_lpr_8(:,2)],'type','samp');
                          
results_all_samp_opt=classify_normals([edge_powers_2(:,1),edge_lpr_2(:,1),...
                              edge_powers_4(:,1),edge_lpr_4(:,1),...
                              edge_powers_8(:,1),edge_lpr_8(:,1)],...
                             [edge_powers_2(:,2),edge_lpr_2(:,2),...
                              edge_powers_4(:,2),edge_lpr_4(:,2),...
                              edge_powers_8(:,2),edge_lpr_8(:,2)],'type','samp',...
                              'dom',results_all.samp_opt_bd,'samp_opt',false);
xlabel('$q_s(${\boldmath$x$}$)$ (2, 4 and 8px)','interpreter','latex');
xlim([-125 150]); set(gca,'ytick',[]);
set(gca,'fontsize',13); box off
