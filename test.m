normals=struct;
normals(1).mu=[0;0];
normals(1).v=[2 .5; .5 1];
normals(1).p=.5;

normals(2).mu=[2;1];
normals(2).v=[1 .5; .5 3];
normals(2).p=.3;

normals(3).mu=[-1;-1];
normals(3).v=[2 0; 0 2];
normals(3).p=.2;

[r,r_sign]=opt_bd_multiclass([0;0],normals);