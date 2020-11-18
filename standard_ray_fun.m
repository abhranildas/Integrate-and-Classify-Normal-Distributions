function f=standard_ray_fun(fun,mu,v,n,r,fun_level)
z=r.*n;
x=num2cell(sqrtm(v)*z+mu,2);
f=fun(x{:})-fun_level;
end