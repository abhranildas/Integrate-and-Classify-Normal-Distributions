normals=struct;
normals(1).mu=[1;0];
normals(1).v=1.01*eye(2);
normals(1).p=.25;

normals(2).mu=[0;1];
normals(2).v=eye(2);
normals(2).p=.25;

normals(3).mu=[-1;0];
normals(3).v=eye(2);
normals(3).p=.25;

normals(4).mu=[0;-1];
normals(4).v=eye(2);
normals(4).p=.25;

dth=pi/1e3;
th=-pi/2:dth:pi/2;
n=[cos(th);sin(th)];
    
[r,r_sign]=opt_bd_multiclass(n,normals);

bd_pts_std=cellfun(@(x,y) x.*y,r,num2cell(n,1),'un',0);
bd_pts_std=horzcat(bd_pts_std{:});
r_sign=horzcat(r_sign{:});
plot(bd_pts_std(1,r_sign==1),bd_pts_std(2,r_sign==1),'.b')
hold on
plot(bd_pts_std(1,r_sign==-1),bd_pts_std(2,r_sign==-1),'.r')