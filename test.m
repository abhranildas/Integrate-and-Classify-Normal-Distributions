%% performance benchmark comparison with Monte Carlo
% d1/d2
mu_1=[2;4]; v_1=[1 1.5; 1.5 3];
mu_2=[5;0]; v_2=[3 0; 0 1];

tic
results=classify_normals([mu_1,v_1],[mu_2,v_2],'method','ray','plotmode',0,'samp_opt',false,'AbsTol',0,'RelTol',0);
toc

axis image; axis([-10 10 -10 10])
bd=results.norm_bd;

linear_bd.q2=zeros(2);
linear_bd.q1=[-.7;1];
linear_bd.q0=0;

tic
results=classify_normals([mu_1,v_1],[mu_2,v_2],'dom',linear_bd,'method','gx2','plotmode',0,'samp_opt',false,'AbsTol',0,'RelTol',1e-10);
toc

n_samp=3e5;
n_iter=10;
pe=nan(n_iter,1);
pe_linear=nan(n_iter,1);

parfor i=1:n_iter
    i
    samp_1=mvnrnd(mu_1,v_1,n_samp);
    samp_2=mvnrnd(mu_2,v_2,n_samp);
    
    pe(i)=1-samp_value(samp_1,samp_2,bd)/(2*n_samp);
    pe_linear(i)=1-samp_value(samp_1,samp_2,linear_bd)/(2*n_samp);
end

% std(pe)/mean(pe)
std(pe_linear)/mean(pe_linear)


tic
n_samp=3e5;
samp_1=mvnrnd(mu_1,v_1,n_samp);
samp_2=mvnrnd(mu_2,v_2,n_samp);
pe=1-samp_value(samp_1,samp_2,linear_bd)/(2*n_samp);
toc

%% performance benchmark comparison with Monte Carlo
% e

n_samp=1e4;
% samp_1_2e=exp(mvnrnd([0 0],eye(2),n_samp));
% samp_2_2e=-exp(mvnrnd([1 1],eye(2),n_samp))+[2 2];

results=classify_normals(samp_1_2e,samp_2_2e,'input_type','samp');

tic
results=classify_normals(samp_1_2e,samp_2_2e,'input_type','samp','method','gx2','samp_opt',false,'plotmode',0,'AbsTol',0,'RelTol',1e-2);
toc

bd=results.norm_bd;

mu_1=mean(samp_1_2e)';
v_1=cov(samp_1_2e);

mu_2=mean(samp_2_2e)';
v_2=cov(samp_2_2e);

n_samp=3e4;
n_iter=40;
pe=nan(n_iter,1);

parfor i=1:n_iter
%     i
    samp_1=mvnrnd(mu_1,v_1,n_samp);
    samp_2=mvnrnd(mu_2,v_2,n_samp);
    
    pe(i)=1-samp_value(samp_1,samp_2,bd)/(2*n_samp);
end

std(pe)/mean(pe)

tic
n_samp=3e4;
samp_1=mvnrnd(mu_1,v_1,n_samp);
samp_2=mvnrnd(mu_2,v_2,n_samp);
pe=1-samp_value(samp_1,samp_2,bd)/(2*n_samp);
toc


%% performance benchmark comparison with Monte Carlo
% g1

mu=[1;1;1;1]; v=diag([1 2 3 4]);

quad.q2=eye(4);
quad.q1=zeros(4,1);
quad.q0=-25;

tic
p=integrate_normal(mu,v,quad,'method','ray','mc_samples',1e6,'plotmode',0)
toc
% xlim([-100 100])



n_samp=5e3;
n_iter=10;
p=nan(n_iter,1);

tic
parfor i=1:n_iter
    i
    p(i)=integrate_normal(mu,v,quad,'method','ray','mc_samples',n_samp,'plotmode',0);
end
toc

std(p)/mean(p)

n_samp=5e3;
tic
p=integrate_normal(mu,v,quad,'method','ray','mc_samples',n_samp,'plotmode',0);
toc

std(p)/mean(p)

tic
parfor i=1:n_iter
    i
    samp=mvnrnd(mu,v,n_samp);
    samp_correct=dot(samp,samp*quad.q2',2) + samp*quad.q1 + quad.q0 > 0;
    p(i)=sum(samp_correct)/n_samp;
end
toc

std(p)/mean(p)

tic
n_samp=1e5;
samp=mvnrnd(mu,v,n_samp);
samp_correct=dot(samp,samp*quad.q2',2) + samp*quad.q1 + quad.q0 > 0;
p=sum(samp_correct)/n_samp;
toc

%% performance benchmark comparison with Monte Carlo
% g2

mu_1=[0;0;0;0];
v_1=[1 0 0 0;
     0 2 -1 0;
     0 -1 3 2;
     0 0 2 4];

mu_2=1.5*[1;1;1;1];
v_2=diag([2 2 2 1]);

tic
results=classify_normals([mu_1,v_1],[mu_2,v_2],'prior_1',.7,'vals',[4 0; 0 1],'plotmode',0,'method','ray','mc_samples',1e5);
toc
% xlim([-20 30])
bd=results.norm_bd;

n_samp=2e3;
n_iter=10;
p=nan(n_iter,1);

tic
parfor i=1:n_iter
    i
    results=classify_normals([mu_1,v_1],[mu_2,v_2],'prior_1',.7,'vals',[4 0; 0 1],'plotmode',0,'method','ray','mc_samples',n_samp);
    p(i)=results.norm_err;
end
toc

std(p)/mean(p)

n_samp=2e3;
tic
results=classify_normals([mu_1,v_1],[mu_2,v_2],'prior_1',.7,'vals',[4 0; 0 1],'plotmode',0,'method','ray','mc_samples',n_samp);
toc

std(p)/mean(p)


n_samp=1e4;
n_iter=10;
pe=nan(n_iter,1);

tic
parfor i=1:n_iter
    i
    samp_1=mvnrnd(mu_1,v_1,n_samp);
    samp_2=mvnrnd(mu_2,v_2,n_samp);
    
    [~,~,samp_1_correct,samp_2_correct]=samp_value(samp_1,samp_2,bd,'vals',[4 0; 0 1]);
    pe(i)= .7*mean(~samp_1_correct) + .3*mean(~samp_2_correct);
end
toc

std(pe)/mean(pe)


tic
n_samp=1e4;
samp_1=mvnrnd(mu_1,v_1,n_samp);
samp_2=mvnrnd(mu_2,v_2,n_samp);
[~,~,samp_1_correct,samp_2_correct]=samp_value(samp_1,samp_2,bd,'vals',[4 0; 0 1]);
pe= .7*mean(~samp_1_correct) + .3*mean(~samp_2_correct);
toc

%% Peformance against increasing d'
mu_1=[0;0;0];

% both normals have the same covariance, so we can check against
% the true d' (Mahalanobis distance).
v=[1 .5 .7;
  .5  2  1 ;
  .7  1  3];

steps=linspace(1,100,40);
d_true=nan(size(steps));
d_gx2=nan(size(steps));
d_ray=nan(size(steps));

% timerec=nan(length(steps),2);
parfor i=1:length(steps)
    i
    mu_2=steps(i)*[1;1;1];
    
%     tic
    results_gx2=classify_normals([mu_1,v],[mu_2,v],'method','gx2','AbsTol',0,'RelTol',0,'plotmode',false);
%     t1=toc;
    
    bd=results_gx2.norm_bd;
    d_true(i)=results_gx2.norm_d_a;
%     d_gx2(i)=results_gx2.norm_d_b;
    
%     tic
%     results_ray=classify_normals([mu_1,v],[mu_2,v],'method','ray','AbsTol',0,'RelTol',0,'plotmode',false);
%     t2=toc;
%     timerec(i,:)=[t1 t2];

%     d_ray(i)=results_ray.norm_d_b;
    
    tic
    samp_1=mvnrnd(mu_1,v,n_samp);
    samp_2=mvnrnd(mu_2,v,n_samp);    
    d_mc(i)=-2*norminv(samp_value(samp_1,samp_2,bd,'vals',[0 1; 1 0])/(2*n_samp));  
    toc
    
end

% for Monte Carlo
steps_mc=[linspace(1,9,20), linspace(9.1,11.5,20)];
d_true_mc=nan(size(steps_mc));
niter=10;
d_mc=nan(niter,length(steps_mc));

nsamp_mc=1e8;

parfor i=1:length(steps_mc)
    mu_2=steps_mc(i)*[1;1;1];
    results_gx2=classify_normals([mu_1,v],[mu_2,v],'method','gx2','plotmode',0);
    bd=results_gx2.norm_bd;
    d_true_mc(i)=results_gx2.norm_d_a;
    
    for iter=1:niter
        [i iter]
        samp_1=mvnrnd(mu_1,v,nsamp_mc);
        samp_2=mvnrnd(mu_2,v,nsamp_mc);
        d_mc(iter,i)=-2*norminv(samp_value(samp_1,samp_2,bd,'vals',[0 1; 1 0])/(2*nsamp_mc));
    end
end

rel_err_gx2=abs(d_gx2-d_true)./d_true;
rel_err_ray=abs(d_ray-d_true)./d_true;
rel_err_mc=abs(d_mc-d_true_mc)./d_true_mc;



figure; hold on
rel_err_gx2(rel_err_gx2==0)=1e-19;
rel_err_ray(rel_err_ray==0)=1e-19;
plot(d_true,rel_err_gx2,'-o','color','g','markersize',3,'markerfacecolor','g')
plot(d_true,rel_err_ray,'-o','color','m','markersize',3,'markerfacecolor','m')
plot(d_true_mc,rel_err_mc,'.','color','k','markersize',10)
xline(-2*norminv(realmin)) % largest computable d', corr. to the smallest possible error representable in double-precision 
yline(eps) % machine epsilon for double precision
xlim([1 80]) 
xlabel 'true d'''
set(gca,'yscale','log','ytick',10.^(-16:5:-1)); ylim([1e-19 1e-1])
ylabel('rel. error')
legend('gen. chi sq.','ray','Monte Carlo')
hold off

% CLEAN THIS AND PUT IT IN DOC

%% perf. benchmark against Monte Carlo: torus
mu=[0;0;0];
v=[1 0 0;
   0 8 4;
   0 4 8];

fun_torus=@(x1,x2,x3) 1.5-(5-(x1.^2+x2.^2).^0.5).^2-x3.^2;
tic
p=integrate_normal(mu,v,fun_torus,'dom_type','fun','fun_span',5,'fun_resol',100,'AbsTol',0,'RelTol',1e-1,'plotmode',0)
toc

n_samp=2e5;
n_iter=10;
p=nan(n_iter,1);

tic
parfor i=1:n_iter
    i
    samp=mvnrnd(mu,v,n_samp);
    samp_cell=num2cell(samp,1);
    p(i)=mean(fun_torus(samp_cell{:})>0);
end
toc

std(p)/mean(p)

n_samp=2e5;
tic
samp=mvnrnd(mu,v,n_samp);
samp_cell=num2cell(samp,1);
p=mean(fun_torus(samp_cell{:})>0);
toc

%% perf. benchmark against Monte Carlo: union of circles
mu=[0;0]; v=[.5 0; 0 1];

left_circle_fun=@(x,y) -(x+1).^2-y.^2+5;
right_circle_fun=@(x,y) -(x-1).^2-y.^2+5;
union_fun=@(x,y) max(left_circle_fun(x,y), right_circle_fun(x,y));

tic
p=integrate_normal(mu,v,union_fun,'dom_type','fun','fun_span',5,'fun_resol',100,'AbsTol',0,'RelTol',1e-2,'plotmode',0)
toc

n_samp=2e2;
n_iter=10;
p=nan(n_iter,1);

tic
parfor i=1:n_iter
    i
    samp=mvnrnd(mu,v,n_samp);
    samp_cell=num2cell(samp,1);
    p(i)=mean(union_fun(samp_cell{:})>0);
end
toc

std(p)/mean(p)

n_samp=2e2;
tic
samp=mvnrnd(mu,v,n_samp);
samp_cell=num2cell(samp,1);
p=mean(union_fun(samp_cell{:})>0);
toc

%%
% mu_1_3d=[0;0;0];
% v_1_3d=eye(3);
% 
% mu_2_3d=[1;1;1];
% v_2_3d=[4    0  -.3  ;
%         0    1  -.5  ;
%       -.3  -.5    2  ];

mu_1_2d=[2;4]; v_1_2d=[1 0; 0 3];
mu_2_2d=[5;0]; v_2_2d=[3 0; 0 .01];

sigma=10.^linspace(-6,2,100);
s_list=10.^linspace(0,3,10);

d_b_1d=nan(length(s_list),length(sigma));
d_a_1d=nan(size(d_b_1d));
d_e_1d=nan(size(d_b_1d));

d_b_2d=nan(size(d_b_1d));
d_a_2d=nan(size(d_b_1d));
d_e_2d=nan(size(d_b_1d));

for i=1:length(s_list)
    s=s_list(i);
    d_a_1d(i,:)=1./(sigma*sqrt((s^2+1)/2));
    d_e_1d(i,:)=2./(sigma*(s+1));
    
    parfor j=1:length(sigma)
        [i j]
        results_1d=classify_normals([0,(s*sigma(j))^2],[1,sigma(j)^2],'plotmode',false);
        d_b_1d(i,j)=results_1d.norm_d_b;
        
        v_1_2d_this=(s*sigma(j))^2*v_1_2d;
        v_2_2d_this=sigma(j)^2*v_2_2d;
        results_2d=classify_normals([mu_1_2d,v_1_2d_this],[mu_2_2d,v_2_2d_this],'plotmode',false);
        d_b_2d(i,j)=results_2d.norm_d_b;
        d_a_2d(i,j)=results_2d.norm_d_a;
        d_e_2d(i,j)=results_2d.norm_d_e;
    end
end

figure;
subplot(1,2,1); hold on; axis([0 60 0 1.15])
xlabel("$d'_b$",'interpreter','latex')
ylabel("$d'_a/d'_b, d'_e/d'_b$",'interpreter','latex')
title('1d')
set(gca,'fontsize',13,'xtick',[0 60])
subplot(1,2,2); hold on; axis([0 60 0 1.15])
title('2d')
set(gca,'fontsize',13,'xtick',[0 60])
cmap=colormap(flipud(winter(length(s_list)))); c=colorbar('ticks',linspace(0,1,4),'ticklabels',[1 10 100 1000]);
set(get(c,'label'),'string','s');

for i=1:size(d_b_1d,1)
    subplot(1,2,1)
    plot(d_b_1d(i,:),d_e_1d(i,:)./d_b_1d(i,:),'color',cmap(i,:))
    plot(d_b_1d(i,:),d_a_1d(i,:)./d_b_1d(i,:),'--','color',cmap(i,:))
    
    subplot(1,2,2)
    plot(d_b_2d(i,:),d_e_2d(i,:)./d_b_2d(i,:),'color',cmap(i,:))
    plot(d_b_2d(i,:),d_a_2d(i,:)./d_b_2d(i,:),'--','color',cmap(i,:))
end

%% Asymptotic properties figure
mu_a=0; s=.05; s=2;
mu_b=1; s_b=.05;
results=classify_normals([mu_a,(s*s)^2],[mu_b,s^2]);
axis([-.5 2.5 0 4.2])
title ''
set(gca,'xtick',[0 1],'ytick',[],'box','off','fontsize',13,'tickdir','out')
%%
colors=colororder;
th=linspace(0,2*pi,100);
z=[cos(th);sin(th)];

mua=[1;2];
mub=[5;-1];
Sa=[2 -.7; -.7 1];
Sb=[.4 -.4; -.4 1.5];
S_avg=(Sa+Sb)/2;
S_rms=sqrtm((Sa^2+Sb^2)/2);

xa=Sa*z+mua;
xb=Sb*z+mub;
x_avg_a=S_avg*z+mua; x_avg_b=S_avg*z+mub;
x_rms_a=S_rms*z+mua; x_rms_b=S_rms*z+mub;

figure; hold on
plot(xa(1,:),xa(2,:),'color',colors(1,:),'linewidth',1)
plot(xb(1,:),xb(2,:),'color',colors(2,:),'linewidth',1)
plot(x_avg_a(1,:),x_avg_a(2,:),'color',colors(5,:),'linewidth',.75)
plot(x_avg_b(1,:),x_avg_b(2,:),'color',colors(5,:),'linewidth',.75)
plot(x_rms_a(1,:),x_rms_a(2,:),'color',colors(3,:),'linewidth',.75)
plot(x_rms_b(1,:),x_rms_b(2,:),'color',colors(3,:),'linewidth',.75)
axis image; axis([-1.5 7 -3 4])
set(gca,'xtick',[],'ytick',[],'fontsize',13,'box','on')

% after whitening and rotation (directly)
S=sqrtm(Sa\Sb^2/Sa);
[R,D]=eig(S);

muaL=[0;0];
mubL=R'/Sa*(mub-mua);
saL=eye(2);
sbL=D;
diag_a=diag(saL); diag_b=diag(sbL);
sL_avmid_a=diag([(diag_a(1)+diag_b(1))/2 diag_a(2)]);
sL_avmid_b=diag([(diag_a(1)+diag_b(1))/2 diag_b(2)]);
sL_avg=(saL+sbL)/2;
sL_rms=sqrtm((saL^2+sbL^2)/2);

xaL=saL*z+muaL; xbL=sbL*z+mubL;
xL_avmid_a=sL_avmid_a*z+muaL; xL_avmid_b=sL_avmid_b*z+mubL;
xL_avg_a=sL_avg*z+muaL; xL_avg_b=sL_avg*z+mubL;
xL_rms_a=sL_rms*z+muaL; xL_rms_b=sL_rms*z+mubL;

figure; hold on
plot(xaL(1,:),xaL(2,:),'color',colors(1,:),'linewidth',1)
plot(xbL(1,:),xbL(2,:),'color',colors(2,:),'linewidth',1)
plot(xL_avmid_a(1,:),xL_avmid_a(2,:),'--','color',colors(5,:),'linewidth',.75)
plot(xL_avmid_b(1,:),xL_avmid_b(2,:),'--','color',colors(5,:),'linewidth',.75)
plot(xL_avg_a(1,:),xL_avg_a(2,:),'color',colors(5,:),'linewidth',.75)
plot(xL_avg_b(1,:),xL_avg_b(2,:),'color',colors(5,:),'linewidth',.75)
plot(xL_rms_a(1,:),xL_rms_a(2,:),'color',colors(3,:),'linewidth',.75)
plot(xL_rms_b(1,:),xL_rms_b(2,:),'color',colors(3,:),'linewidth',.75)
axis image; axis([-3 1.5 -2 4])
set(gca,'xtick',[],'ytick',[],'fontsize',13,'box','on')

xaM=S_avg\xa;
xbM=S_avg\xb;
x_avg_aM=S_avg\x_avg_a;
x_avg_bM=S_avg\x_avg_b;
figure; hold on
plot(xaM(1,:),xaM(2,:),'color',colors(1,:),'linewidth',1)
plot(xbM(1,:),xbM(2,:),'color',colors(2,:),'linewidth',1)
plot(x_avg_aM(1,:),x_avg_aM(2,:),'color',colors(5,:),'linewidth',1)
plot(x_avg_bM(1,:),x_avg_bM(2,:),'color',colors(5,:),'linewidth',1)

%% yes/no and 2I task
ma=2.3; va=1;
mb=0; vb=1;

results_yn=classify_normals([ma,va],[mb,vb]);

m_ab=[ma; mb];
v_ab=diag([va vb]);

m_ba=[mb; ma];
v_ba=diag([vb va]);

results_2I=classify_normals([m_ab,v_ab],[m_ba,v_ba]);

axis image; axis([-5 5 -5 5])

%% multi-interval task
pc=[];
for n=2:9 % # of intervals
    
    qm.q2=1/vb-1/va;
    qm.q1=2*(ma/va-mb/vb);
    qm.q0=mb^2/vb-ma^2/va;
    
    [w_a,k_a,lambda_a,m_a,s_a]=gx2_params_norm_quad(ma,va,qm);
%     [fa,xa]=gx2pdf('full',w_a,k_a,lambda_a,m_a,s_a);
    
    [w_b,k_b,lambda_b,m_b,s_b]=gx2_params_norm_quad(mb,vb,qm);
%     [fb,xb]=gx2pdf('full',w_b,k_b,lambda_b,m_b,s_b);
    
    % figure; hold on
    % plot(xa,fa)
    % plot(xb,fb)
    % axis([-1 15 1e-2 1e2])
    % set(gca,'yscale','log')
    
%     q=linspace(-10,20,1e3);
%     dq=q(2)-q(1);
    [mu,v]=gx2stat(w_a,k_a,lambda_a,m_a,s_a);
    fun=@(q)gx2pdf(q,w_a,k_a,lambda_a,m_a,s_a).*...
        gx2cdf(q,w_b,k_b,lambda_b,m_b,s_b).^(n-1);
%     figure; plot(q,f)
%     pc=integral(fun,-10,20)
    pc=[pc; [n integral(fun,mu-5*sqrt(v),mu+5*sqrt(v))]];
end
plot(pc(:,1),pc(:,2),'-o')
xlabel 'number of intervals'
ylabel 'accuracy'

%% ROC curves figure
ma=0; va=1;
mb=1; vb=4;

ma_4d=[0;0;0;0];
va_4d=eye(4);

mb_4d=.2*[1;1;-6;0];
vb_4d=diag([4 4 2 1]);

sc=linspace(-4,10,100);
lr=linspace(-1,10,100);
lr_4d=linspace(-15,10,100);

results=classify_normals([ma,va],[mb,vb],'plotmode',0);
hf=results.norm_errmat(:,2)*2;
d_a=results.norm_d_a;
quad_lr=results.norm_bd;

results_4d=classify_normals([ma_4d,va_4d],[mb_4d,vb_4d],'plotmode','fun_prob');
hf_4d=results_4d.norm_errmat(:,2)*2;
quad_lr_4d=results_4d.norm_bd;

quad_sc.q2=0;
quad_sc.q1=-1;

ROC_sc=nan(length(sc),2);
ROC_lr=nan(length(lr),2);
ROC_lr_4d=nan(length(lr_4d),2);

for i=1:length(sc)
    quad_sc.q0=sc(i);
    results=classify_normals([ma,va],[mb,vb],'dom',quad_sc,'plotmode',0);
    ROC_sc(i,:)=results.norm_errmat(:,2)'*2;
    
    quad_lr.q0=lr(i);
    results=classify_normals([ma,va],[mb,vb],'dom',quad_lr,'plotmode',0);
    ROC_lr(i,:)=results.norm_errmat(:,2)'*2;
    
    quad_lr_4d.q0=lr_4d(i);
    results=classify_normals([ma_4d,va_4d],[mb_4d,vb_4d],'method','gx2','dom',quad_lr_4d,'plotmode',0);
    ROC_lr_4d(i,:)=results.norm_errmat(:,2)'*2;
end

figure; hold on
plot(ROC_sc(:,1),ROC_sc(:,2))
plot(ROC_lr(:,1),ROC_lr(:,2))
plot(hf(1),hf(2),'o')
% plot(ROC_lr_4d(:,1),ROC_lr_4d(:,2))
% plot(hf_4d(1),hf_4d(2),'o')
axis image; axis([0 1 0 1]);
r=refline(1,0); set(r,'color','k')
xlabel('p(f)'); ylabel('p(h)');

%AUC_sc=abs(trapz(ROC_sc(:,1),ROC_sc(:,2)));
%[d_a sqrt(2)*norminv(AUC_sc)]
AUC_lr=abs(trapz(ROC_lr(:,1),ROC_lr(:,2)));

%% get d'_b from ROC curve
n_samp=20;

ROC_sc=nan(length(sc),2);
for i=1:length(sc)
    samp_a=normrnd(ma,sqrt(va),[n_samp 1]);
    samp_b=normrnd(mb,sqrt(vb),[n_samp 1]);
    quad_sc.q0=sc(i);
    results=classify_normals(samp_a,samp_b,'input_type','samp','dom',quad_sc,'plotmode',0,'samp_opt',false);
    ROC_sc(i,:)=results.norm_errmat(:,2)'*2;    
end

ROC_z=norminv(ROC_sc);
p=polyfit(ROC_z(:,1),ROC_z(:,2),1);
plot(ROC_z(:,1),ROC_z(:,2),'o')
refline(p(1),p(2))
axis image

m_hat=p(2)/p(1)
s_hat=1/p(1)

results=classify_normals([0,1],[m_hat,s_hat^2])
%%
n_samp=1e6;

ma=[2;0;0;0];
va=eye(4);

% mb=[1;1;1;1];
vb=diag([1 3 1 5]);
t_dof=5;

xa=mvnrnd(ma,va,n_samp);
xb=mvtrnd(vb,t_dof,n_samp);
% xb=mvnrnd(mb,vb,n_samp);

la=log(mvtpdf(xa,vb,t_dof))-log(mvnpdf(xa,ma',va));
lb=log(mvtpdf(xb,vb,t_dof))-log(mvnpdf(xb,ma',va));

cleanidx=(~isnan(la)&~isnan(lb)&~isinf(la)&~isinf(lb));
la=la(cleanidx);
lb=lb(cleanidx);

histogram(la,'normalization','pdf','edgecolor','none'); hold on
histogram(lb,'normalization','pdf','edgecolor','none')
xlim([-5 20])

ldiff=lb-la;
% figure; histogram(ldiff,'normalization','pdf','edgecolor','none')
% xlim([-5 20])
p_2I=mean(ldiff>0)

% results=classify_normals(la,lb,'input_type','samp','samp_opt',false);

l=[linspace(-10,30,200) linspace(40,700,10)];
ROC=nan(length(l),2);

quad.q2=0;
quad.q1=-1;
for i=1:length(l)
    quad.q0=l(i);
    results=classify_normals(la,lb,'input_type','samp','samp_opt',false,'dom',quad,'plotmode',0);
    ROC(i,:)=results.samp_errmat(:,1)'/length(la);
end

figure; hold on
plot(ROC(:,2),ROC(:,1))
axis image; axis([0 1 0 1]);
r=refline(1,0); set(r,'color','k')
xlabel('p(f)'); ylabel('p(h)');

AUC=abs(trapz(ROC(:,2),ROC(:,1)))

%% log normal
mu=1; v=.25;

fun=@(x) sin(exp(x));
integrate_normal(mu,v,fun,'dom_type','fun','fun_span',10,'fun_resol',500)

x=linspace(-2,2,100);
f=norm_fun_pdf(x,mu,v,fun,'fun_span',10,'fun_resol',500);
plot(x,f); ylabel 'cdf'