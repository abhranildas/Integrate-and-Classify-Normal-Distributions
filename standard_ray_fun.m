function f=standard_ray_fun(fun,mu,v,n,r)
% dim=length(mu);
z=r.*n;
%x=;
% if dim==1
%     f=func(x);
% elseif dim==2
%     f=func(x(1,:),x(2,:));
% elseif dim==3
%     f=func(x(1,:),x(2,:),x(3,:));
% end
x=num2cell(sqrtm(v)*z+mu,2);
f=fun(x{:});
end