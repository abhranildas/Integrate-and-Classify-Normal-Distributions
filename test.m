mu=[0;0];
v=eye(2);

cheb_reg_x=@(x,y) y;
reg_1=@(n,orig)ray_scan(cheb_reg_x,n,orig);

cheb_reg_2=@(x,y) x;
reg_2=@(n,orig)ray_scan(cheb_reg_2,n,orig);

reg_intersect=@(n,orig) combine_regs({reg_1,reg_2},'or',n,orig);
integrate_normal(mu,v,@(n) reg_intersect(n,mu),'reg_type','ray_scan','n_rays',1e3);
axis image
xlim([-10 10])
ylim([-10 10])
%%
mu=[0;0;0];
v=eye(3);
orig=[1;1;1];

reg_x=@(n,orig)ray_scan(@(x,y,z) x,n,orig);

reg_y=@(n,orig)ray_scan(@(x,y,z) y,n,orig);

reg_z=@(n,orig)ray_scan(@(x,y,z) z,n,orig);

reg_vol_min=@(n,orig)ray_scan(@(x,y,z) x*y*z-1,n,orig);

reg_vol_max=@(n,orig)ray_scan(@(x,y,z) 5-x*y*z,n,orig);

reg_surf_min=@(n,orig)ray_scan(@(x,y,z) 2*(x*y+y*z+z*x)-10,n,orig);

reg_surf_max=@(n,orig)ray_scan(@(x,y,z) 12-2*(x*y+y*z+z*x),n,orig);

reg_intersect=@(n,orig) combine_regs({reg_vol_min,reg_vol_max,reg_surf_min,reg_surf_max},'and',n,orig);

plot_boundary(reg_intersect,3,'n_rays',1e3,'orig',orig)

xlim([0 5])
ylim([0 5])
zlim([0 5])

cheb_reg_2=@(x,y) x;
reg_2=@(n,orig)ray_scan(cheb_reg_2,n,orig);

reg_intersect=@(n,orig) combine_regs({reg_1,reg_2},'or',n,orig);

integrate_normal(mu,v,@(n) reg_1(n,mu),'reg_type','ray_scan','n_rays',1e2);
%integrate_normal(mu,v,@(n) reg_intersect(n,mu),'reg_type','ray_scan','n_rays',1e3);
axis image
xlim([-10 10])
ylim([-10 10])

%%
mu=[0;0];
v=eye(2);
orig=[1;1];

reg_y_pos=@(n,orig)ray_scan(@(x,y) y,n,orig);

reg_damped_osc=@(n,orig)ray_scan(@(x,y) abs(exp(-abs(x))*sin(x))-y-.2,n,orig);

reg_intersect=@(n,orig) combine_regs({reg_y_pos,reg_damped_osc},'and',n,orig);

plot_boundary(reg_damped_osc,2,'n_rays',1e1,'orig',orig)

xlim([0 5])
ylim([0 5])
zlim([0 5])

cheb_reg_2=@(x,y) x;
reg_2=@(n,orig)ray_scan(cheb_reg_2,n,orig);

reg_intersect=@(n,orig) combine_regs({reg_1,reg_2},'or',n,orig);

integrate_normal(mu,v,@(n) reg_1(n,mu),'reg_type','ray_scan','n_rays',1e2);
%integrate_normal(mu,v,@(n) reg_intersect(n,mu),'reg_type','ray_scan','n_rays',1e3);
axis image
xlim([-10 10])
ylim([-10 10])

%%
d = [-3 10 -3 3];
f = chebfun2('y.*cos(y.^2+x)-.1',d);
g = chebfun2('cos(x.^2/2).*sin(y.^2)-.1',d);
% Plot zero contours of f & g
froot=roots(f);
plot(roots(f)), hold on, plot(roots(g))
% Plot their common roots
r = roots(f, g, 'resultant');
plot(r(:,1), r(:,2), 'o')

e=1;
cheb.xy;
f=x.*y.^2/2-e; % define the region using a cheb function
froot=roots(f);

f = chebfun2(@(x,y) x.*y.^2/2-1,[inf inf -inf inf]);
froot=roots(f);
plot(froot)

%%
mu=[0;0;0];
v=eye(3);

quad.a2=-eye(3);
quad.a1=[0;0;0];
quad.a0=100;

[p,pc]=int_norm_grid(mu,v,@(n) ray_scan(quad,'quad',n,mu))


%%
lambda=[-10 -5 7];
m=[4 20 3];
delta=[200 3 10];

gx2cdf_imhof(0,lambda,m,delta,'tail')

x=linspace(-5e3,-1e3,1e2);
p=arrayfun(@(x) gx2cdf_imhof(x,lambda,m,delta),x);
p_tail=arrayfun(@(x) gx2cdf_imhof(x,lambda,m,delta,'tail'),x);
figure(1)
plot(x,p,'-');
hold on
plot(x,p_tail,'o');
hold off

%% symbolic integration
lambda=[1 5 2];
m=[2 3 7];
delta=[1 2 3];
x=1;

syms u

F=vpaintegral(@(u) integrand(u,x,lambda',m',delta'),u,0,inf)

gx2cdf_imhof(1,[1 5 2],[1 2 3],[2 3 7])

function f=integrand(u,x,lambda,m,delta)
% lambda=[1 5 2]';
% m=[1 2 3]';
% delta=[2 3 7]';
% x=1;
theta=sum(m.*atan(lambda*u)+(delta.*(lambda*u))./(1+lambda.^2*u.^2),1)/2-u*x/2;
rho=prod(((1+lambda.^2*u.^2).^(m/4)).*exp(((lambda.^2*u.^2).*delta)./(2*(1+lambda.^2*u.^2))),1);
f=sin(theta)./(u.*rho);
% f=u^2;
end

%     function f=imhof_integrand(u,x,lambda,m,delta)
%         theta=sum(m.*atan(lambda*u)+(delta.*(lambda*u))./(1+lambda.^2*u.^2),1)/2-u*x/2;
%         rho=prod(((1+lambda.^2*u.^2).^(m/4)).*exp(((lambda.^2*u.^2).*delta)./(2*(1+lambda.^2*u.^2))),1);
%         f=sin(theta)./(u.*rho);
%     end