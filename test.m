%%
mu_1=[4; 5]; v_1=[2 1; 1 1];
mu_2=mu_1; v_2=3*[2 -1; -1 1];

results_ray=classify_normals([mu_1,v_1],[mu_2,v_2],'plotmode','fun_prob')

%%
mu_1=0; v_1=1;
mu_2=2.5; v_2=1.5;

fun=@(x) cos(x.^2);
results=classify_normals([mu_1,v_1],[mu_2,v_2],'dom',fun,'dom_type','fun','plotmode','fun_prob')

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

%% conv test
% x = [5 6 8 2 5];
% y = [6 -1 3 5];
% x1 = [x zeros(1,3)]
% y1 = [y zeros(1,4)]
% c1 = ifft(fft(x1).*fft(y1))
% c2 = conv(x,y)

% w=[1 -5];
% k=[1 2];
% lambda=[2 3];
% m=0;
% s=0;

% % current gx2pdf
% [mu,v]=gx2stat(w,k,lambda,m,s);
% xg=linspace(mu-5*sqrt(v),mu+5*sqrt(v),1e3);
% f=gx2pdf(xg,w,k,lambda,m,s);
% figure; plot(xg,f)
% xlim([min(xg) max(xg)])

% % ncx2 pdf's
% [mu1,v1]=gx2stat(w(1),k(1),lambda(1),0,0);
% [mu2,v2]=gx2stat(w(2),k(2),lambda(2),0,0);
% span=max([abs(mu1)+5*sqrt(v1) abs(mu2)+5*sqrt(v2)]);
% x=linspace(-span,span,2e3);
% dx=diff(x); dx=dx(1);
% 
% f1=gx2pdf(x,w(1),k(1),lambda(1),0,0);
% f2=gx2pdf(x,w(2),k(2),lambda(2),0,0);
% figure; hold on
% plot(x,f1); plot(x,f2);

% c1 = ifftshift(ifft(fft(f1).*fft(f2)))*dx;
% hold on; plot(x,c1); 
% 
% c2=conv(f1,f2)*dx;
% xc2=linspace(-2*span+dx,2*span-dx,length(c2));
% figure; plot(xc2,c2)
% xlim([min(xg) max(xg)])

% w=[-5 2 3];
% k=[1 2 100];
% lambda=[2 3 7];
% m=-50;
% s=4;
% 
% [fplot,x]=gx2pdf('plot',w,k,lambda,m,s);
% plot(x,fplot);
% 
% x=[-500 200 275 300 623];
% f=gx2pdf(x,w,k,lambda,m,s);
% hold on
% plot(x,f,'o')


mu=[1;1;1;1]; v=diag([1 2 3 4]);

quad.q2=eye(4);
quad.q1=zeros(4,1);
quad.q0=-25;
[w,k,lambda,m,s]=gx2_params_norm_quad(mu,v,quad);

[f_conv,x]=gx2pdf('plot',w,k,lambda,m,s);
plot(x,f_conv);
f_diff=gx2pdf(x,w,k,lambda,m,s,'method','diff');
hold on;
plot(x,f_diff)

%% 2I task
mu_n=2; sig_n=1;
mu_s=3; sig_s=2;

mu_a=[mu_s; mu_n];
v_a=diag([sig_s^2 sig_n^2]);

mu_b=[mu_n; mu_s];
v_b=diag([sig_n^2 sig_s^2]);

results=classify_normals([mu_a,v_a],[mu_b,v_b]);
axis image; axis([0 5.5 0 5.5])
r=refline(1,0);
set(r,'color','g')

d1=(mu_s-mu_n)*sig_s/(sig_s^2-sig_n^2);
d2=(mu_s-mu_n)*sig_n/(sig_s^2-sig_n^2);
gx2cdf(0,[sig_s^2 -sig_n^2],[1 1], [d1^2 d2^2],0,0)

[w,k,lambda,s,m]=gx2_params_norm_quad(mu_a,v_a,results.norm_bd)
gx2cdf(0,w,k,lambda,s,m)