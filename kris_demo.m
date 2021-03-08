% circle at (1,1) with radius 0.1
dom.q2=eye(2);
dom.q1=-[2;2];
dom.q0=1.9;

% scan it with rays
n=[1 1;0 1];
[init_sign,x]=ray_scan(dom,n)

% normal parameters
mu=[-1;-1];
v=eye(2);

% now scan domain standardized wrt normal
[init_sign,x]=ray_scan(dom,n,'mu',mu,'v',v)

% see the whole picture
integrate_normal(mu,v,dom)

%% ray-scan demo
mu=[1;7];
v=[1 .5; .5 2];

integrate_normal(mu,v,@line_ray_scan,'dom_type','ray_scan')
axis([-10 10 -10 10])

% line y>1 in quad form, i.e. 0.x^2+0.y^2+0.xy+0.x+1.y-1>0,
% i.e. [0 0; 0 0]*[x;y] + [0;1]*[x;y] -1 >0
line_quad.q2=zeros(2);
line_quad.q1=[0;1];
line_quad.q0=-1;

figure;
integrate_normal(mu,v,line_quad,'dom_type','quad')
axis([-10 10 -10 10])

% function form
line_fun=@(x,y) y-1;

figure;
integrate_normal(mu,v,line_fun,'dom_type','fun')
axis([-10 10 -10 10])

mu_2=[-1;-2];
v_2=[1 -.5; -.5 2];
n_samp=1e3;
samp_1=mvnrnd(mu,v,n_samp);
samp_2=mvnrnd(mu_2,v_2,n_samp);
results=classify_normals(samp_1,samp_2,'input_type','samp','dom',@line_ray_scan,'dom_type','ray_scan');
axis([-10 10 -10 10])