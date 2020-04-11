mu=[0;0];
v=eye(2);

cheb_reg_x=@(x,y) y;
reg_1=@(n,orig)ray_scan(cheb_reg_x,n,orig);

cheb_reg_2=@(x,y) x;
reg_2=@(n,orig)ray_scan(cheb_reg_2,n,orig);

reg_intersect=@(n,orig) combine_regs({reg_1,reg_2},'or',n,orig);
integrate_normal(mu,v,@(n) reg_intersect(n,mu),'reg_type','ray_scan');
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

%% y>sin(x)
mu=[0;0];
v=eye(2);

reg_cheb=@(x,y) y-exp(-x.^2).*sin(10*x);
reg_rayscan=@(n,orig)ray_scan(reg_cheb,'cheb',n,orig);
plot_ray_scan_bd(reg_rayscan,2,'n_rays',1e3)

integrate_normal(mu,v,reg_cheb,'reg_type','cheb','n_bd_pts',1e3);
xlim([-10 10])
ylim([-15 15])
