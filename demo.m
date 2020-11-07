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

% integrate in a quadratic region x^2-x-1<0, i.e. f(x)= -1*x^2 + 1*x + 1 >0
reg_quad.q2=-1;
reg_quad.q1=1;
reg_quad.q0=1;
integrate_normal(mu,v,reg_quad);

% most plots can be zoomed and panned.

% integrate in a region defined by a non-quadratic f(x)>0
fun=@(x) cos(x.^2);
figure;
integrate_normal(mu,v,fun,'reg_type','fun','fun_span',5,'fun_resol',500);

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

% integrate in a quadratic region (x+y)^2 > x+1, i.e.
% q(x,y)= [x y]*[1 1; 1 1]*[x;y] + [-1 0]*[x;y] -1 > 0
reg_quad.q2=[1 1; 1 1];
reg_quad.q1=[-1;0];
reg_quad.q0=-1;

% compare two integration algorithms
figure; integrate_normal(mu,v,reg_quad); % ray method
figure; integrate_normal(mu,v,reg_quad,'method','gx2'); % gx2 method

%% 2D, classify two normals, one inside the other
mu_1=[4; 5];
v_1=[2 1; 1 1];

mu_2=mu_1;
v_2=3*[2 -1; -1 1];

% compare two integration algorithms
results_ray=classify_normals([mu_1,v_1],[mu_2,v_2])
results_gx2=classify_normals([mu_1,v_1],[mu_2,v_2],'method','gx2')

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
plot_boundary(linear_bd,2,'reg_type','quad','mu',[-4;4],'plot_type','line','line_color',[1 0 1])
axis image; xlim([-10 10]); ylim([-10 10])
set(gca,'fontsize',13); box off

results_linear=classify_normals([mu_1,v_1],[mu_2,v_2],'reg',linear_bd,'bPlot',false)

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

results_samp_custom=classify_normals(samp_1,samp_2,'type','samp','reg',custom_bd)
axis image; xlim([-10 10]); ylim([-10 10])

%% PAPER 2D, classify non-normal samples
n_samp=1e3;
samp_1=exp(mvnrnd([0 0],eye(2),n_samp));
samp_2=-exp(mvnrnd([1 1],eye(2),n_samp));

results=classify_normals(samp_1,samp_2,'type','samp')
axis image; xlim([-20 15]); ylim([-20 15])
set(gca,'fontsize',13); box off

%% Inversion/union/intersection of integration/classification regions
mu=[0;0];
v=[.5 0; 0 1];

% circle_left.q2=-eye(2);
% circle_left.q1=[-2;0];
% circle_left.q0=4;
% 
% circle_right=circle_left;
% circle_right.q1=-circle_left.q1;
% 
% circle_left_rayscan=@(n,mu,v)ray_scan(circle_left,n,'mu',mu,'v',v);
% circle_right_rayscan=@(n,mu,v)ray_scan(circle_right,n,'mu',mu,'v',v);

% circle_union_rayscan=@(n,mu,v) combine_regs({circle_left_rayscan,circle_right_rayscan},'or',n,'mu',mu,'v',v);
% circle_intersection_rayscan=@(n,mu,v) combine_regs({circle_left_rayscan,circle_right_rayscan},'and',n,'mu',mu,'v',v);

% integrate_normal(mu,v,circle_union_rayscan,'reg_type','ray_scan','AbsTol',0,'RelTol',1e-10);
% axis image; xlim([-4 4]); ylim([-4 4])

circle_left=@(x,y) -(x+1).^2-y.^2+5;
circle_right=@(x,y) -(x-1).^2-y.^2+5;

circle_union=@(x,y) max(circle_left(x,y), circle_right(x,y));
integrate_normal(mu,v,circle_union,'reg_type','fun','fun_span',5);
axis image; xlim([-4 4]); ylim([-4 4])

% figure
% integrate_normal(mu,v,circle_intersection_rayscan,'reg_type','ray_scan');
% axis image; xlim([-4 4]); ylim([-4 4])

circle_intersection=@(x,y) min(circle_left(x,y), circle_right(x,y));

figure
integrate_normal(mu,v,circle_intersection,'reg_type','fun','fun_span',5);
axis image; xlim([-4 4]); ylim([-4 4])

% circle_right_invert_rayscan=@(n,mu,v) invert_reg(circle_right_rayscan,n,'mu',mu,'v',v);

% crescent_rayscan=@(n,mu,v) combine_regs({circle_left_rayscan,circle_right_invert_rayscan},'or',n,'mu',mu,'v',v);
% 
% figure
% integrate_normal(mu,v,crescent_rayscan,'reg_type','ray_scan');
% axis image; xlim([-4 4]); ylim([-4 4])

crescent=@(x,y) max(circle_left(x,y), -circle_right(x,y));

figure
integrate_normal(mu,v,crescent,'reg_type','fun','fun_span',5);
axis image; xlim([-4 4]); ylim([-4 4])

% classify normals using this region
mu_2=[2.2;0];
v_2=[.5 0; 0 .25];
classify_normals([mu,v],[mu_2,v_2],'reg',crescent,'reg_type','fun','fun_span',5);
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
    
    results_gx2=classify_normals([mu_1,v],[mu_2,v],'method','gx2','AbsTol',0,'RelTol',0,'bPlot',false);
    d_true(i)=results_gx2.norm_maha_dprime;
    d_gx2(i)=results_gx2.norm_dprime;
    
    results_ray=classify_normals([mu_1,v],[mu_2,v],'method','ray','AbsTol',0,'RelTol',0,'bPlot',false);
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

%% PAPER Integrate in a toroidal region defined by implicit function f(x)>0

mu=[0;0;0];
v=[1 0 0;
   0 8 4;
   0 4 8];

% plot the error ellipsoid of the normal
plot_normal(mu,v);

fun_torus=@(x1,x2,x3) 1.5-(5-(x1.^2+x2.^2).^0.5).^2-x3.^2;
figure;
integrate_normal(mu,v,fun_torus,'reg_type','fun','fun_span',3,'fun_resol',10,'RelTol',1e-1);
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
reg_quad.q2=eye(4);
reg_quad.q1=zeros(4,1);
reg_quad.q0=-25;

p=integrate_normal(mu,v,reg_quad)
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
tic
results=classify_normals([mu_1,v_1],[mu_2,v_2],'prior_1',.7,'vals',[4 0; 0 1],'method','ray','mc_samples',5e2,'bPlot',false)
t=toc
ylim([0 .12])

% now classify using samples
n_samp=1e3;
results=classify_normals(mvnrnd(mu_1',v_1,n_samp),mvnrnd(mu_2',v_2,n_samp),'type','samp','prior_1',.7,'vals',[4 0; 0 1])
set(gca,'fontsize',13); box off

%% PAPER 4D Integrate in a custom region, using Monte Carlo
mu=zeros(4,1);
v=eye(4);

fun=@(x1,x2,x3,x4) 1-abs(x1)-abs(x2)-abs(x3)-abs(x4);
mc_samples=round(10.^linspace(1,4,20)); % # of Monte Carlo samples
n_repeat=5;
plist=nan(length(mc_samples),n_repeat);

for k=1:n_repeat
    for i=1:length(mc_samples)
        [k i]
        p=integrate_normal(mu,v,fun,'reg_type','fun','fun_span',3,'fun_resol',10,'mc_samples',mc_samples(i))
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
integrate_normal(mu,v,f1,'reg_type','fun');
axis image; xlim([-10 20]); ylim([-50 10]);

% integrate f2 above 0
figure;
integrate_normal(mu,v,f2,'reg_type','fun');
axis image; xlim([-10 20]); ylim([-50 10]);

% prob. that cos(x) and cos(y) are both >0
% i.e. integrate vector function f=[f1,f2] in the domain f1>0 and f2>0, i.e. min(f1,f2)>0.
f_domain1=@(x,y) min(f1(x,y),f2(x,y));

figure;
integrate_normal(mu,v,f_domain1,'reg_type','fun','fun_span',10);
axis image; xlim([-10 20]); ylim([-50 10]);
set(gca,'fontsize',13); box off

% prob. that cos(x) + cos(y) > 0.5
% i.e. integrate the vector function in the implicit domain f1+f2-0.5 > 0
f_domain2=@(x,y) f1(x,y)+f2(x,y)-.5;

figure;
integrate_normal(mu,v,f_domain2,'reg_type','fun','fun_span',10);
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
plot_boundary(@(n,mu,v) opt_class_multi(n,normals,3,'mu',mu,'v',v),2,'mu',normals(3).mu,'reg_type','ray_scan')
title 'boundary of normal 3'

% classify
results=classify_normals_multi(normals)
axis image

q12=opt_class_quad([normals(1).mu normals(1).v],[normals(2).mu normals(2).v]);
q13=opt_class_quad([normals(1).mu normals(1).v],[normals(3).mu normals(3).v]);
f12=quad2fun(q12,1);
f13=quad2fun(q13,1);
f=@(x,y) min(f12(x,y),f13(x,y));
plot_boundary(f,2,'reg_type','fun','plot_type','line')

% now classify using samples from these normals
samples=struct;
for i=1:4
    samples(i).sample=mvnrnd(normals(i).mu,normals(i).v,1e4);
end
results_samp=classify_normals_multi(samples,'type','samp')
axis image

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
normals=struct;
normals(1).mu=[0;0;-1;-1]; normals(1).v=diag([1 2 3 4]);
normals(2).mu=[1;1;0;0]; normals(2).v=eye(4)/4;
normals(3).mu=[2;2;1;1]; normals(3).v=eye(4);
normals(4).mu=[3;3;3;3]; normals(4).v=eye(4);

priors=[.3 .3 .35 .05];
results=classify_normals_multi(normals,'priors',priors,'mc_samples',1e3,'plotmode',[1;1;1;1])

% now classify using samples from these normals
samples=struct;
for i=1:4
    samples(i).sample=mvnrnd(normals(i).mu,normals(i).v,2e2);
end
results_samp=classify_normals_multi(samples,'type','samp','priors',priors,'mc_samples',1e3,'plotmode',[1;1;1;1])
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

results_dv_2=classify_normals(results_joint_2.samp_opt_dv{1},results_joint_2.samp_opt_dv{2},'type','samp','reg',@(x) x,'reg_type','fun','samp_opt',false);
xlim([-50 100]); ylim([0 .04]); set(gca,'ytick',[]);
set(gca,'fontsize',13); box off
xlabel('$q_s(${\boldmath$x$}$)$','interpreter','latex');

results_joint_4=classify_normals([edge_powers_4(:,1),edge_lpr_4(:,1)],[edge_powers_4(:,2),edge_lpr_4(:,2)],'type','samp','bPlot',false);
results_joint_8=classify_normals([edge_powers_8(:,1),edge_lpr_8(:,1)],[edge_powers_8(:,2),edge_lpr_8(:,2)],'type','samp','bPlot',false);

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
                              'reg',results_all.samp_opt_bd,'samp_opt',false);
xlabel('$q_s(${\boldmath$x$}$)$ (2, 4 and 8px)','interpreter','latex');
xlim([-125 150]); set(gca,'ytick',[]);
set(gca,'fontsize',13); box off
