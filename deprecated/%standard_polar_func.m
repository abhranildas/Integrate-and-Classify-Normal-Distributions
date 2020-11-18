function f=standard_polar_func(func,mu,v,r,az,el)
dim=length(mu);

if dim==1
    f=func(sqrt(v)*r+mu);
elseif dim==2
    [z1,z2]=pol2cart(az,r);
    x=sqrtm(v)*[z1;z2]+mu;
    f=func(x(1,:),x(2,:));
elseif dim==3
    [z1,z2,z3]=sph2cart(az,el,r);
    x=sqrtm(v)*[z1(:)';z2(:)';z3(:)']+mu;
    f=reshape(func(x(1,:),x(2,:),x(3,:)),size(z1));
end
end