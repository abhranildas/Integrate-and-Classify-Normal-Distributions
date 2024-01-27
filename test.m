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

results=classify_normals([mu_1,v_1],[mu_2,v_2],'eff',0.1)

%% test efficiency scalar: 2d
mu_1=[2;4]; v_1=[1 1.5; 1.5 3];
mu_2=[5;0]; v_2=[3 0; 0 1];
efflist=linspace(1,0,10);
for i=1:length(efflist)
    results=classify_normals([mu_1,v_1],[mu_2,v_2],'eff',efflist(i))
    axis([0 10 -3 7])
    drawnow
    pause
end

%% test efficiency scalar: 2d data
%% make d' scaling animation

mu_1=[0;4]; v_1=[1 1.5; 1.5 3];
mu_2=[5;-1]; v_2=[3 0; 0 2];

n_samp=5e3;
samp_1=mvnrnd(mu_1,v_1,n_samp);
samp_2=mvnrnd(mu_2,v_2,n_samp);

efflist=linspace(1,0,200);

frames = {};
for i=1:length(efflist)
    i
    results=classify_normals(samp_1,samp_2,'input_type','samp','method','gx2','eff',efflist(i),'samp_opt',false);
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

%% animation: scale d' of non-normal samples

mu_1=[-10;5]; v_1=[10 -2; -2 2]*.5;
mu_2=[5;5]; v_2=[1 1.5; 1.5 3]*.1;

n_samp=5e3;
samp_1=mvnrnd(mu_1,v_1,n_samp);
samp_2=mvnrnd(mu_2,v_2,n_samp);

samp_1(:,2)=samp_1(:,2).^2/1.5;
samp_1=samp_1+[-5 10];
samp_2=samp_2.^6/1e3;

efflist=linspace(1,0,200);

frames = {};
for i=1:length(efflist)
    i
    results=classify_normals(samp_1,samp_2,'input_type','samp','method','gx2','eff',efflist(i),'samp_opt',false);
    axis image; axis([-30 60 -20 100])
    box on
    set(gca,'xtick',[],'ytick',[])
%     title(sprintf("$d' = %.1f$",results.norm_d_b),'interpreter','latex');
    text(30,-10,sprintf("$d' = %.1f$",results.norm_d_b),'interpreter','latex','fontsize',30);
%     set(gca,'fontsize',13)
    frame = getframe(gca);
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

