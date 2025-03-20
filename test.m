samp_1=[normrnd(0,2,5e2,1); normrnd(4,1,1e3,1)];
samp_2=[normrnd(2,1,1e3,1); normrnd(6,2,5e2,1)];
results=classify_normals(samp_1,samp_2,'input_type','samp','samp_opt',false,'vals',[2 0; 0 1]);
axis([-6 12 0 .15])

%%
mu_1=[2;2]; v_1=[1 1.5; 1.5 3];
mu_2=[3;0]; v_2=[3 0; 0 1];

n_samp=1e3;
samp_1=mvnrnd(mu_1,v_1,n_samp);
samp_2=mvnrnd(mu_2,v_2,n_samp);

results=classify_normals(samp_1,samp_2,'input_type','samp','samp_opt',false,'vals',[4 0; 0 1])

%%
absent=importdata('target_absent.txt',',',1);
present=importdata('target_present.txt',',',1);

results=classify_normals(present.data,absent.data,'input_type','samp','d_con',true)

%%
normals=struct;
normals(1).mu=[1;0]; normals(1).v=.1*eye(2);
normals(2).mu=[0;1]; normals(2).v=.1*eye(2);
normals(3).mu=[-1;0]; normals(3).v=.1*eye(2);
normals(4).mu=[0;-1]; normals(4).v=.1*eye(2);

samples=struct;
for i=1:4
    samples(i).sample=mvnrnd(normals(i).mu,normals(i).v,1e4);
end
vals=diag([1 2 3 4]);
results=classify_normals_multi(samples,'input_type','samp','vals',vals)

%% test efficiency scalar: 1d
mu_1=0; v_1=1;
mu_2=2.5; v_2=1.5;

results=classify_normals([mu_1,v_1],[mu_2,v_2],'d_scale',0.5)

%% test efficiency scalar: 2d
mu_1=[2;4]; v_1=[1 1.5; 1.5 3];
mu_2=[5;0]; v_2=[3 0; 0 1];
d_scale_list=linspace(1,0,10);
for i=1:length(d_scale_list)
    results=classify_normals([mu_1,v_1],[mu_2,v_2],'d_scale',d_scale_list(i))
    axis([0 10 -3 7])
    drawnow
    pause
end

%% test efficiency scalar: 2d data
%% d' scaling animation

mu_1=[0;4]; v_1=[1 1.5; 1.5 3];
mu_2=[5;-1]; v_2=[3 0; 0 2];

n_samp=5e3;
samp_1=mvnrnd(mu_1,v_1,n_samp);
samp_2=mvnrnd(mu_2,v_2,n_samp);

d_scale_list=linspace(1,0,200);

frames = {};
for i=1:length(d_scale_list)
    i
    results=classify_normals(samp_1,samp_2,'input_type','samp','method','gx2','eff',d_scale_list(i),'samp_opt',false);
    axis image; axis([-5 10 -5 10])
    box on
    set(gca,'xtick',[],'ytick',[])
    title(sprintf("$d' = %.1f$",results.norm_d_b),'interpreter','latex');
    frame = getframe(gcf);
    frames{end+1} = frame;
    close
end

% Create the GIF
gifFileName = 'output.gif';
for i = [1:numel(frames)-1 numel(frames)-1:-1:2]
    im = frame2im(frames{i});
    [A,map] = rgb2ind(im,256);
    if i == 1
        imwrite(A,map,gifFileName,'gif','LoopCount',Inf,'DelayTime',0.03);
    else
        imwrite(A,map,gifFileName,'gif','WriteMode','append','DelayTime',0.03);
    end
end

%% scale d' of non-normal samples

mu_1=[-10;5]; v_1=[10 -2; -2 2]*.5;
mu_2=[5;5]; v_2=[1 1.5; 1.5 3]*.1;

n_samp=5e3;
samp_1=mvnrnd(mu_1,v_1,n_samp);
samp_2=mvnrnd(mu_2,v_2,n_samp);

samp_1(:,2)=samp_1(:,2).^2/1.5;
samp_1=samp_1+[-5 10];
samp_2=samp_2.^6/1e3;

results=classify_normals(samp_1,samp_2,'input_type','samp','method','gx2','d_scale',0.5,'d_scale_type','squeeze_dist','samp_opt',false);
axis image; axis([-30 50 -20 100])

%% scale d' of non-normal samples: animation

mu_1=[-10;5]; v_1=[10 -2; -2 2]*.5;
mu_2=[5;5]; v_2=[1 1.5; 1.5 3]*.1;

n_samp=5e3;
samp_1=mvnrnd(mu_1,v_1,n_samp);
samp_2=mvnrnd(mu_2,v_2,n_samp);

samp_1(:,2)=samp_1(:,2).^2/1.5;
samp_1=samp_1+[-5 10];
samp_2=samp_2.^6/1e3;

d_scale_list=linspace(1,0,100);

frames = {};
for i=1:length(d_scale_list)
    i
    results=classify_normals(samp_1,samp_2,'input_type','samp','method','gx2','d_scale',d_scale_list(i),'samp_opt',false);
    axis image; axis([-30 60 -20 100])
    box on
    set(gca,'xtick',[],'ytick',[])
    %     title(sprintf("$d' = %.1f$",results.norm_d_b),'interpreter','latex');
    title ''
    text(20,-10,sprintf("$d' = %.0f$",results.norm_d_b),'interpreter','latex','fontsize',20);
    %     set(gca,'fontsize',13)
    frame = getframe(gca);
    frames{end+1} = frame;
    % pause
    close
end

% Create the GIF
gifFileName = 'output.gif';
for i = [1:numel(frames)-1 numel(frames)-1:-1:2]
    im = frame2im(frames{i});
    [A,map] = rgb2ind(im,256);
    if i == 1
        imwrite(A,map,gifFileName,'gif','LoopCount',Inf,'DelayTime',0.03);
    else
        imwrite(A,map,gifFileName,'gif','WriteMode','append','DelayTime',0.03);
    end
end

%% scale d' of non-normal samples: figure for paper
%% warping the distributions

mu_1=[-10;5]; v_1=[10 -2; -2 2]*.3;
mu_2=[6;5]; v_2=[1 1.5; 1.5 3]*.02;

n_samp=5e3;
samp_1=mvnrnd(mu_1,v_1,n_samp);
samp_2=mvnrnd(mu_2,v_2,n_samp);

samp_1(:,2)=samp_1(:,2).^2/1.5;
samp_1=samp_1+[-5 10];
samp_2=samp_2.^6/1e3;

d_scale_list=linspace(1,0.05,4);
d_list=nan(size(d_scale_list));
for i=1:length(d_scale_list)
    results=classify_normals(samp_1,samp_2,'input_type','samp','d_scale',d_scale_list(i),'samp_opt',false);
    d_list(i)=results.samp_d_b;
    set(gca,'xtick',[],'ytick',[])
end
axis image; axis([-25 70 0 50])

%% warping the dv

mu_1=[-10;5]; v_1=[10 -2; -2 2]*.3;
mu_2=[6;5]; v_2=[1 1.5; 1.5 3]*.02;

n_samp=5e3;
samp_1=mvnrnd(mu_1,v_1,n_samp);
samp_2=mvnrnd(mu_2,v_2,n_samp);

samp_1(:,2)=samp_1(:,2).^2/1.5;
samp_1=samp_1+[-5 10];
samp_2=samp_2.^6/1e3;

d_scale_list=linspace(0.5,0,4);
d_list=nan(size(d_scale_list));
for i=1:length(d_scale_list)
    % figure
    results=classify_normals(samp_1,samp_2,'input_type','samp','d_scale',d_scale_list(i),'d_scale_type','squeeze_dv','samp_opt',false,'plotmode','fun_prob');
    % subplot(4,1,i); hold on
    % histogram(results.samp_dv{1},'edgecolor','none','Normalization','pdf')
    % histogram(results.samp_dv{2},'edgecolor','none','Normalization','pdf')
    % xline(0)
    xlim([-1200 500])
    % set(gca,'xtick',[-1200 0 500],'ytick',[])
    d_list(i)=results.samp_d_b;
end

%% d'-scaling each dimension
mu_1=[2;4]; v_1=[1 1.5; 1.5 3]; S_1=sqrtm(v_1);
mu_2=[5;0]; v_2=[3 -1.5; -1.5 1]; S_2=sqrtm(v_2);

results=classify_normals([mu_1,v_1],[mu_2,v_2]);
axis([0 10 -3 7])

% first scale x. need to scale both variance and covariance
mu_2(1)=mu_1(1);
scale_x=sqrt(v_1(1,1)/v_2(1,1));
v_2(1,1)=v_2(1,1)*scale_x^2;
v_2(1,2)=v_2(1,2)*scale_x;
v_2(2,1)=v_2(2,1)*scale_x;
results=classify_normals([mu_1,v_1],[mu_2,v_2]);
axis([0 10 -3 7])

% then scale y
mu_2(2)=mu_1(2);
scale_y=sqrt(v_1(2,2)/v_2(2,2));
v_2(2,2)=v_2(2,2)*scale_y^2;
v_2(1,2)=v_2(1,2)*scale_y;
v_2(2,1)=v_2(2,1)*scale_y;
results=classify_normals([mu_1,v_1],[mu_2,v_2]);
axis([0 10 -3 7])

%% pulling bayes dv's together
mu_1=[2;4]; v_1=[1 1.5; 1.5 3];
mu_2=[5;0]; v_2=[3 -1.5; -1.5 1];

n_samp=5e3;
samp_1=mvnrnd(mu_1,v_1,n_samp);
samp_2=mvnrnd(mu_2,v_2,n_samp);

results=classify_normals(samp_1,samp_2,'input_type','samp','samp_opt',false,'d_scale',0);

%% test flipping

mu=[-1; -1]; v=[1 0.5; 0.5 2];

quad.q2=[1 1; 1 1];
quad.q1=[-1;0];
quad.q0=-1;

p=integrate_normal(mu,v,quad,'method','gx2') % gx2 method

figure
[p,~,bd_pts]=integrate_normal(mu,v,quad,'method','ray','add_bd_pts',true)
% [p,~,bd_pts]=int_norm_ray(mu,v,quad)

% mu=[0;0;0];
% v=[ 1    0   .1;
%     0    2    1;
%    .1    1    4];

L=[.5 -1]; U=[1 3];
p=mvncdf(L,U,mu',v)

figure
p=integrate_normal(mu,v,[L;U],'dom_type','rect','add_bd_pts',true)

%% non-orthogonal contributions to d'
% mu_1=[2;3];
% mu_2=[7;9];
% v=[1 -1.3; -1.3 4];

mu_1=[0;0];
mu_2=[9;9];
v=[1 .8; .8 1];

S=sqrtm(v);

samp_1=mvnrnd(mu_1,v,1e3);
samp_2=mvnrnd(mu_2,v,1e3);

classify_normals(samp_1,samp_2,'input_type','samp','samp_opt',false);
hold on
% axis normal

% plot line between means
plot([mu_1(1), mu_2(1)], [mu_1(2), mu_2(2)], '-k')

% plot rectangular grid
s=sqrt(diag(v));
x_min=mu_1(1)-s(1);
x_max=mu_2(1)+s(1);
y_min=mu_1(2)-s(2);
y_max=mu_2(2)+s(2);

x_grid = x_min:s(1):x_max; % X-coordinates for vertical lines
y_grid = y_min:s(2):y_max; % Y-coordinates for horizontal lines

% Generate endpoints for vertical lines
X_vertical = [x_grid; x_grid]; % Each column is a vertical line
Y_vertical = [y_min * ones(size(x_grid)); y_max * ones(size(x_grid))];

% Generate endpoints for horizontal lines
Y_horizontal = [y_grid; y_grid]; % Each column is a horizontal line
X_horizontal = [x_min * ones(size(y_grid)); x_max * ones(size(y_grid))];

% Plot all vertical and horizontal lines at once
% plot(X_vertical, Y_vertical, 'k-', X_horizontal, Y_horizontal, 'k-');

% plot quarter unit circles
th=linspace(0,2*pi,100);
circ=[sin(th);cos(th)];
circ_x=sqrt(diag(v)).*circ+mu_1;
circ_z=S*circ+mu_1;

plot(circ_z(1,:),circ_z(2,:),'-r')
plot(circ_x(1,:),circ_x(2,:),'-k')

% axis([-1 10 -3 15])
title ''
% set(gca,'xtick',[],'ytick',[],'fontsize',13)

hold on

% plot basis vectors of z
quiver(mu_1(1), mu_1(2), S(1,1), S(2,1), 0, '-r','showarrowhead',0);
quiver(mu_1(1), mu_1(2), S(1,2), S(2,2), 0, '-r','showarrowhead',0);

% whitened space:
samp_1_w=samp_1/S;
samp_2_w=samp_2/S;

mu_1_w=S\mu_1;
mu_2_w=S\mu_2;
v_w=eye(2);

results=classify_normals(samp_1_w,samp_2_w,'input_type','samp','samp_opt',false);

T=inv(S);
hold on

% d' vector
dprime_w=mu_2_w-mu_1_w;
quiver(mu_1_w(1), mu_1_w(2), dprime_w(1), dprime_w(2), 0, 'k','showarrowhead',0);

% projections of d' along original axes
d_1=dot(dprime_w,T(:,1))*T(:,1)/norm(T(:,1))^2;
d_2=dot(dprime_w,T(:,2))*T(:,2)/norm(T(:,2))^2;

% individual dprimes
d=mu_2-mu_1;
dprimes_ind=d./sqrt(diag(v))
dprimes_ind(2)/dprimes_ind(1);

% vector contribution lengths
% delta=d.*(vecnorm(T,2,1))'
delta=d.*sqrt(diag(inv(v)))
delta(2)/delta(1);
c=d;

% plot parallelogram
quiver(mu_1_w(1), mu_1_w(2), c(1)*T(1,1), c(1)*T(2,1), 0, 'k','showarrowhead',0);
quiver(mu_1_w(1), mu_1_w(2), c(2)*T(1,2), c(2)*T(2,2), 0, 'k','showarrowhead',0);

quiver(mu_1_w(1)+c(1)*T(1,1), mu_1_w(2)+c(1)*T(2,1), c(2)*T(1,2), c(2)*T(2,2), 0, 'k','showarrowhead',0);
quiver(mu_1_w(1)+c(2)*T(1,2), mu_1_w(2)+c(2)*T(2,2), c(1)*T(1,1), c(1)*T(2,1), 0, 'k','showarrowhead',0);


% plot original basis
quiver(mu_1_w(1), mu_1_w(2), T(1,1), T(2,1), 0, 'b','showarrowhead',0);
quiver(mu_1_w(1), mu_1_w(2), T(1,2), T(2,2), 0, 'b','showarrowhead',0);

% angles between d' vector and original axes
acosd(sum(dprime_w.*T)./(norm(dprime_w)*vecnorm(T)))

title ''
% axis([0 14.5 -1 10.5])
% set(gca,'xtick',[],'ytick',[],'fontsize',13)

%% d_con_vec
mu_1=[0;0;0];
mu_2=[1;1;1];
v=[1 0  0;
  0  1   0;
  0  0   1];
results=classify_normals([mu_1 v],[mu_2 v]);
S=sqrtm(v);
dprime_w=S\(mu_2-mu_1);
T=inv(S);

sum(dprime_w.*T)./(norm(dprime_w)*vecnorm(T))
