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

% you can zoom and pan all plots.

% integrate in a region defined by a non-quadratic f(x)>0
f=@(x) cos(x.^2);
figure;
integrate_normal(mu,v,f,'reg_type','fun','fun_span',5,'fun_resol',500);

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

%% 1D, classify using samples (priors are taken to be prop. to sample sizes)
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

integrate_normal(mu,v,reg_quad);
%xlim([-5 5])

%% 2D, classify two normals, one inside the other
mu_1=[4; 5];
v_1=[2 1; 1 1];

mu_2=mu_1;
v_2=3*[2 -1; -1 1];

results=classify_normals([mu_1,v_1],[mu_2,v_2])

%% PAPER 2D, classify with custom boundaries and from samples

mu_1=[2;4];
v_1=[1 1.5; 1.5 3];

mu_2=[5;0];
v_2=[3 0; 0 1];

results=classify_normals([mu_1,v_1],[mu_2,v_2])

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
samp_2=exp(mvnrnd([1 1],eye(2),n_samp));

results=classify_normals(samp_1,samp_2,'type','samp')
axis image; xlim([-10 20]); ylim([-5 20])
set(gca,'fontsize',13); box off

%% PAPER Inversion/union/intersection of integration/classification regions
mu=[0;0];
v=[.5 0; 0 1];

circle_left.q2=-eye(2);
circle_left.q1=[-2;0];
circle_left.q0=4;

circle_right=circle_left;
circle_right.q1=-circle_left.q1;

circle_left_rayscan=@(n,mu,v)ray_scan(circle_left,n,'mu',mu,'v',v);
circle_right_rayscan=@(n,mu,v)ray_scan(circle_right,n,'mu',mu,'v',v);

circle_union_rayscan=@(n,mu,v) combine_regs({circle_left_rayscan,circle_right_rayscan},'or',n,'mu',mu,'v',v);
circle_intersection_rayscan=@(n,mu,v) combine_regs({circle_left_rayscan,circle_right_rayscan},'and',n,'mu',mu,'v',v);

integrate_normal(mu,v,circle_union_rayscan,'reg_type','ray_scan');
axis image; xlim([-4 4]); ylim([-4 4])

figure
integrate_normal(mu,v,circle_intersection_rayscan,'reg_type','ray_scan');
axis image; xlim([-4 4]); ylim([-4 4])

circle_left_invert_rayscan=@(n,mu,v) invert_reg(circle_left_rayscan,n,'mu',mu,'v',v);

crescent_rayscan=@(n,mu,v) combine_regs({circle_left_invert_rayscan,circle_right_rayscan},'and',n,'mu',mu,'v',v);

figure
integrate_normal(mu,v,crescent_rayscan,'reg_type','ray_scan');
axis image; xlim([-4 4]); ylim([-4 4])

% classify normals using this region
mu_2=[2.2;0];
v_2=[.5 0; 0 .25];
classify_normals([mu,v],[mu_2,v_2],'reg',crescent_rayscan,'reg_type','ray_scan');
axis image; xlim([-2 4]); ylim([-3 3])
set(gca,'fontsize',13); box off

%% 3D, classify
mu_1=[0;0;0];
v_1=eye(3);

mu_2=[2;1;1];
v_2=2*eye(3);

results=classify_normals([mu_1,v_1],[mu_2,v_2])

% plot the boundary points used for integration
%hold on
%plot3(results.norm_bd_pts(1,:),results.norm_bd_pts(2,:),results.norm_bd_pts(3,:),'.','markersize',4)

%% High-accuracy estimation of tiny errors (large d')
format long
dprime_true=75

mu_1=[0;0;0];
v_1=eye(3);

mu_2=dprime_true*[1;0;0];
v_2=(1+1e-12)*eye(3);

results=classify_normals([mu_1,v_1],[mu_2,v_2],'AbsTol',0,'RelTol',1e-3)
axis image; xlim([-5 80]); ylim([-2 2]); zlim([-2 2]); view(-53,25)
dprime_computed=results.norm_dprime
format

%% PAPER Integrate in a toroidal region defined by implicit function f(x)>0

mu=[0;0;0];
v=[1 0 0;
   0 8 4;
   0 4 8];

torus=@(x1,x2,x3) 1.5-(5-(x1.^2+x2.^2).^0.5).^2-x3.^2;
integrate_normal(mu,v,torus,'reg_type','fun','fun_span',5,'fun_resol',10);
axis image; xlim([-7 7]); ylim([-7 7]); zlim([-7 7]);
set(gca,'xtick',linspace(-7,7,5)); set(gca,'ytick',linspace(-7,7,5)); set(gca,'ztick',linspace(-7,7,5))
set(gca,'fontsize',13)

%% PAPER 4D, integrate
mu=[1;1;1;1];
v=diag([1 2 3 4]);

% integrate outside the sphere of radius 5, i.e.
% x1^2+x2^2+x3^2+x4^2-5^2>0, i.e. x'*eye(4)*x + zeros(4,1)'*x -25 >0
reg_quad.q2=eye(4);
reg_quad.q1=zeros(4,1);
reg_quad.q0=-25;

integrate_normal(mu,v,reg_quad);
xlim([-40 40]); ylim([0 .06]);
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
ylim([0 .12])

% now classify using samples
n_samp=1e3;
results=classify_normals(mvnrnd(mu_1',v_1,n_samp),mvnrnd(mu_2',v_2,n_samp),'type','samp','prior_1',.7,'vals',[4 0; 0 1])
set(gca,'fontsize',13); box off

%% REMOVE THIS
% Integrate non-quadratic function f of a normal,
% equivalent to integrating normal in the non-quadratic region f>0

mu_1=[2;3];
v_1=[1 -.5; -.5 2];
e=12;

e_fun=@(x,y) e-x.*y.^2/2; % define f

integrate_normal(mu_1,v_1,e_fun,'reg_type','fun');

% classify normals and samples wrt this region
mu_2=[3;4];
v_2=[2 -1; -1 4];
n_samp=1e3;
samp_1=mvnrnd(mu_1,v_1,n_samp);
samp_2=mvnrnd(mu_2,v_2,n_samp);
results=classify_normals(samp_1,samp_2,'prior_1',.8,'type','samp','reg',e_fun,'reg_type','fun');

% integrate vector-valued function of a normal
% integrate [f,g]=[e_func,p_func] where both are +ve
p=7;
p_quad.q2=(eye(2)-1)/2;
p_quad.q1=[0;0];
p_quad.q0=p;

e_rayscan=@(n,mu,v)ray_scan(e_fun,n,'mu',mu,'v',v,'reg_type','fun','fun_span',15);
p_rayscan=@(n,mu,v)ray_scan(p_quad,n,'mu',mu,'v',v);
ep_rayscan=@(n,mu,v) combine_regs({e_rayscan,p_rayscan},'and',n,'mu',mu,'v',v);

figure; hold on
plot_boundary(e_fun,2,'reg_type','fun','plot_color','b')
plot_boundary(p_quad,2,'reg_type','quad','plot_color','r')
xlim([-5 15])
ylim([-15 15])

% classify normals and samples using this region
results=classify_normals([mu_1,v_1],[mu_2,v_2],'prior_1',.8,'reg',ep_rayscan,'reg_type','ray_scan');
xlim([-5 20]); ylim([-15 15])
results=classify_normals(samp_1,samp_2,'prior_1',.8,'type','samp','reg',ep_rayscan,'reg_type','ray_scan');
xlim([-5 20]); ylim([-15 15])

%% PAPER Integrate vector-valued function of a normal

mu=[6;-19];
v=[1 -.7; -.7 2];

f1=@(x,y) cos(x);
f2=@(x,y) cos(y);

% integrate f1 above 0
figure;
integrate_normal(mu,v,f1,'reg_type','fun');

% integrate f2 above 0
figure;
integrate_normal(mu,v,f2,'reg_type','fun');

% first convert the regions to ray-scan format
% f1r=@(n,mu,v)ray_scan(f1,n,'mu',mu,'v',v,'reg_type','func','func_span',15);
% f2r=@(n,mu,v)ray_scan(f2,n,'mu',mu,'v',v,'reg_type','func','func_span',15);

% you can combine ray-scan regions, here using 'and', i.e. intersection.
% f_domain1=@(n,mu,v) combine_regs({f1r,f2r},'and',n,'mu',mu,'v',v);

% now integrate in the combined region
% figure;
% integrate_normal(mu,v,f_domain1,'reg_type','ray_scan');

% integrate vector function f=[f1,f2] in the domain f1>0 and f2>0, i.e. min(f1,f2)>0.
f_domain1=@(x,y) min(f1(x,y),f2(x,y));

figure;
integrate_normal(mu,v,f_domain1,'reg_type','fun','fun_span',10);
axis image; xlim([-10 20]); ylim([-50 10]);
set(gca,'fontsize',13); box off

% integrate the vector function in an implicit domain such as f1+f2 > 0.2
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

% define means and vcovs of the normals
mus=[[1;0],[0;1],[-1;0],[0;-1]];
vs=cat(3,2*eye(2),eye(2),eye(2),eye(2));

% define struct of all normals
normals=struct;
for i=1:4
    normals(i).mu=mus(:,i); normals(i).v=vs(:,:,i);
end

% plot the multi-class boundary of normal 3
plot_boundary(@(n,mu,v) opt_class_multi(n,mus,vs,3,'mu',mu,'v',v),2,'mu',mus(:,3),'reg_type','ray_scan')
title 'boundary of normal 3'

% classify
results=classify_normals_multi(normals)
axis image

% now classify using samples from these normals
dists2=struct;
for i=1:4
    dists2(i).sample=mvnrnd(normals(i).mu,normals(i).v,1e2);
end
results_samp=classify_normals_multi(dists2,'type','samp')
axis image

%% PAPER Classifying 7 normals, 2D
% define means and vcovs of the normals
mus=[[2;0],[0;1],[-1;0],[0;-1],[-2;2],[2;-3],[-2;-2.5]];
vs=cat(3,...
    [2 1; 1 2],...
    [.5 0; 0 1],...
    [1 .3; .3 1],...
    [.5 -.5; -.5 1],...
    .5*[1 1; 1 5],...
    .3*eye(2),...
    [1 0; 0 .1]);

% define struct of all normals
normals=struct;
for i=1:7
    normals(i).mu=mus(:,i); normals(i).v=vs(:,:,i);
end

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

%% PAPER Actual vision research data: occluding target detection

dataArray = textscan(fopen('absent.txt','r'), '%*q%f%f%f%[^\n\r]', 'Delimiter', ',', 'HeaderLines' ,1);
absent = [dataArray{1:end-1}];

dataArray = textscan(fopen('present.txt','r'), '%*q%f%f%f%[^\n\r]', 'Delimiter', ',', 'HeaderLines' ,1);
present = [dataArray{1:end-1}];

results=classify_normals(absent,present,'type','samp')
axis normal
xlim([-1 2]); ylim([0 .5]); zlim([-400 400]); view(28,18);
set(gca,'ytick',0:.2:.4);
xlabel('edge'); ylabel('luminance'); zlabel('pattern');
set(gca,'fontsize',13); box off

%% PAPER Actual vision research data: camouflage detection

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
