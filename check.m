%% no quad term

% 1D
mu=1;
v=2;

% 2x>3 -> 2x-3>0
reg_quad.q2=0;
reg_quad.q1=2;
reg_quad.q0=-3;

p_ray=integrate_normal(mu,v,reg_quad)
p_gx2=int_norm_quad_gx2(mu,v,reg_quad)

% 2D
mu=[0;0];
v=eye(2);

% x+y>1 -> x+y-1>0
reg_quad.q2=[0 0; 0 0];
reg_quad.q1=[1;1];
reg_quad.q0=-1;

p_ray=integrate_normal(mu,v,reg_quad)
xlim([-5 5])
p_gx2=int_norm_quad_gx2(mu,v,reg_quad)

% 3D
mu=[1;2;3];
v=2*eye(3);

% x+2y+3z>1
reg_quad.q2=zeros(3);
reg_quad.q1=[1;2;3];
reg_quad.q0=-1;

p_ray=integrate_normal(mu,v,reg_quad)
xlim([-5 5])
p_gx2=int_norm_quad_gx2(mu,v,reg_quad)

%% parabola
% 2D
mu=[0;0];
v=eye(2);

% y^2-x>0
reg_quad.q2=[0 0; 0 1];
reg_quad.q1=[-1;0];
reg_quad.q0=0;

% (y-1)^2>(x-1) -> y^2-x-2y+2>0
reg_quad.q2=[0 0; 0 1];
reg_quad.q1=[-1;-2];
reg_quad.q0=2;

p_ray=integrate_normal(mu,v,reg_quad)
xlim([-5 5])
p_gx2=int_norm_quad_gx2(mu,v,reg_quad)

% 3D
mu=[1;2;3];
v=[1 0 0; 0 2 0; 0 0 3];

% elliptic paraboloid: 2x^2+3y^2-z>0
reg_quad.q2=[2 0 0; 0 3 0; 0 0 0];
reg_quad.q1=[0;0;-1];
reg_quad.q0=0;

p_ray=integrate_normal(mu,v,reg_quad)
xlim([-5 5]); ylim([-5 5]); zlim([-5 5])
p_gx2=int_norm_quad_gx2(mu,v,reg_quad)

% hyperbolic paraboloid: 2x^2-3y^2-z>0
reg_quad.q2=[2 0 0; 0 -3 0; 0 0 0];
reg_quad.q1=[0;0;-1];
reg_quad.q0=0;

p_ray=integrate_normal(mu,v,reg_quad,'AbsTol',1e-10,'RelTol',1e-3)
xlim([-5 5]); ylim([-5 5]); zlim([-5 5])
p_gx2=int_norm_quad_gx2(mu,v,reg_quad,'AbsTol',1e-10,'RelTol',1e-3)

%% point: all/no space

% 1D
mu=2;
v=3;

% (x-1)^2>0
reg_quad.q2=1;
reg_quad.q1=-2;
reg_quad.q0=1;

p_ray=integrate_normal(mu,v,reg_quad)
xlim([-5 5])
p_gx2=int_norm_quad_gx2(mu,v,reg_quad)

% 2D
mu=[0;0];
v=eye(2);

% (x-1)^2+(y-1)^2>0 -> x^2 + y^2 - 2xy - 2x - 2y + 2 >0
reg_quad.q2=[2 -1; -1 2];
reg_quad.q1=[-2;-2];
reg_quad.q0=2;

integrate_normal(mu,v,reg_quad);
xlim([-5 5])
p=int_norm_quad_gx2(mu,v,reg_quad)

% 3D
mu=[2;3;5];
v=2*eye(3);

% x^2+y^2+z^2<0
reg_quad.q2=-eye(3);
reg_quad.q1=[0;0;0];
reg_quad.q0=0;

p_ray=integrate_normal(mu,v,reg_quad)
xlim([-5 5])
p_gx2=int_norm_quad_gx2(mu,v,reg_quad)

% 4D
mu=[2;3;5;7];
v=2*eye(4);

% x^2+y^2+z^2<0
reg_quad.q2=-eye(4);
reg_quad.q1=[0;0;0;0];
reg_quad.q0=0;

p_ray=integrate_normal(mu,v,reg_quad)
xlim([-5 5])
p_gx2=int_norm_quad_gx2(mu,v,reg_quad)

%% point: all/no space
mu=[0;0];
v=eye(2);

% (x-1)^2+(y-1)^2+(x-y)^2>0 -> 2x^2 + 2y^2 - 2xy - 2x - 2y + 2 >0
reg_quad.a2=[2 -1; -1 2];
reg_quad.a1=[-2;-2];
reg_quad.a0=2;

% (x-y+1)^2 + (3x+2y-2)^2>0 -> x^2+y^2-2xy+2x-2y+1 +
% 9x^2+4y^2+12xy-12x-8y+4>0 -> 10x^2+5y^2+10xy-10x-10y+5>0
reg_quad.a2=[10 5; 5 5];
reg_quad.a1=[-10;-10];
reg_quad.a0=5;

integrate_normal(mu,v,reg_quad)

p=int_norm_quad_gx2(mu,v,reg_quad)

%% point: all/no space
mu=[0;0;0];
v=eye(3);

% (x-1)^2+(y-1)^2+(z-1)^2>0 -> x^2 + y^2 +z^2 - 2x - 2y - 2z + 3 >0
reg_quad.a2=eye(3);
reg_quad.a1=-2*[1;1;1];
reg_quad.a0=3;

integrate_normal(mu,v,reg_quad);

p=int_norm_quad_gx2(mu,v,reg_quad)

%% intersecting lines
% 2D
mu=[.5;.1];
v=eye(2);

% x1^2-x2^2 > 0
reg_quad.q2=[1 0; 0 -1];
reg_quad.q1=[0;0];
reg_quad.q0=0;

p_ray=integrate_normal(mu,v,reg_quad)
xlim([-5 5])
p_gx2=int_norm_quad_gx2(mu,v,reg_quad)

% 3D
mu=[.5;.1;2];
v=1.5*eye(3);

% (x+y+z)(x+y-z)>0 -> x^2+y^2+2xy-z^2>0
reg_quad.q2=[1 1 0; 1 1 0; 0 0 -1];
reg_quad.q1=[0;0;0];
reg_quad.q0=0;

p_ray=integrate_normal(mu,v,reg_quad)
xlim([-5 5])
p_gx2=int_norm_quad_gx2(mu,v,reg_quad)

%% parallel lines

% 2D
mu=[-2;0];
v=eye(2);

% y=x and y=x+1 -> (y-x)(y-x-1)>0 -> x^2+y^2-2xy+x-y>0
reg_quad.q2=[1 -1; -1 1];
reg_quad.q1=[1;-1];
reg_quad.q0=0;

p_ray=integrate_normal(mu,v,reg_quad)
xlim([-5 5])
p=int_norm_quad_gx2(mu,v,reg_quad)

% 3D
mu=[0;1;0];
v=2*eye(3);

% (x+2y+3z)(x+2y+3z+1) ->
reg_quad.q2=[1 2 3 ;
             2 4 6 ;
             3 6 9];
reg_quad.q1=[1;2;3];
reg_quad.q0=0;

p_ray=integrate_normal(mu,v,reg_quad)
xlim([-5 5])
p=int_norm_quad_gx2(mu,v,reg_quad)

%% coincident/nonexistent lines
% 2D
mu=[0;0];
v=eye(2);

% y=2x+3 -> (y-2x-3)^2 > 0 -> 4x^2+y^2-4xy+12x-6y+9 > 0
reg_quad.a2=[4 -2; -2 1];
reg_quad.a1=[12;-6];
reg_quad.a0=9;

integrate_normal(mu,v,reg_quad);
xlim([-5 5])

p=int_norm_quad_gx2(mu,v,reg_quad)

% 3D
mu=[0;0;0];
v=eye(3);

% (x+3y-2z)^2
reg_quad.q2=[1 3 -2 ;
             3 9 -6 ;
             -2 -6 4];
reg_quad.q1=[0;0;0];
reg_quad.q0=0;

p_ray=integrate_normal(mu,v,reg_quad)
xlim([-5 5])
p_gx2=int_norm_quad_gx2(mu,v,reg_quad)