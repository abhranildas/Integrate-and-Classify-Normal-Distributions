normals=struct;
normals(1).mu=[1;0];
normals(1).v=eye(2);
normals(1).prior=.25;

normals(2).mu=[0;1];
normals(2).v=eye(2);
normals(2).prior=.25;

normals(3).mu=[-1;0];
normals(3).v=eye(2);
normals(3).prior=.25;

normals(4).mu=[0;-1];
normals(4).v=eye(2);
normals(4).prior=.25;

normals=normals(circshift(1:4,-1));

dth=pi/1e3;
th=-pi/2:dth:pi/2;
n=[cos(th);sin(th)];
    
[r,r_sign]=opt_bd_multi(n,normals,'dist_wrt',[normals(2).mu,normals(2).v]);

bd_pts_std=cellfun(@(x,y) x.*y,r,num2cell(n,1),'un',0);
bd_pts_std=horzcat(bd_pts_std{:});
r_sign=horzcat(r_sign{:});
plot(bd_pts_std(1,r_sign==1),bd_pts_std(2,r_sign==1),'.b')
hold on
plot(bd_pts_std(1,r_sign==-1),bd_pts_std(2,r_sign==-1),'.r')

% [p,pc,bd_pts]=integrate_normal(normals(1).mu,normals(1).v,'bd_fn',@(n) opt_bd_multi(n,normals,'id',2));
results=classify_normals_multi(normals);

%%
mu_a=[-1;0];
v_a=eye(2);

mu_b=[1;0];
v_b=2*eye(2);

mu_wrt=[0;2];
v_wrt=eye(2);

dth=pi/1e3;
th=-pi/2:dth:pi/2;
n=[cos(th);sin(th)];
    
[r,r_sign]=opt_bd(n,[mu_a,v_a],[mu_b,v_b],'dist_wrt',[mu_wrt,v_wrt]);

bd_pts_std=cellfun(@(x,y) x.*y,r,num2cell(n,1),'un',0);
bd_pts_std=horzcat(bd_pts_std{:});
r_sign=horzcat(r_sign{:});
plot(bd_pts_std(1,r_sign==1),bd_pts_std(2,r_sign==1),'.b')
hold on
plot(bd_pts_std(1,r_sign==-1),bd_pts_std(2,r_sign==-1),'.r')


%%
mu_a=[0;0;0];
mu_b=[60;0;0];
v=eye(3);

% force grid method
bd_fn_a=@(n) opt_bd(n,[mu_a,v],[mu_b,v]);
bd_fn_b=@(n) opt_bd(n,[mu_b,v],[mu_a,v]);

results_grid=classify_normals([mu_a,v],[mu_b,v],'custom_bd_fns',{bd_fn_a,bd_fn_b});

