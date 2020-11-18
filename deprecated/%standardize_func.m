function func_s=standardize_func(func,mu,v)
dim=length(mu);
    function fs=func_s_array(z)
        x=num2cell(sqrtm(v)*z+mu);
        fs=func(x{:});
    end

if dim==1
    func_s=@(z) sqrtm(v).*z+mu;
elseif dim==2
    func_s=@(z1,z2) arrayfun(@(z1,z2) func_s_array([z1;z2]),z1,z2);
elseif dim==3
    func_s=@(z1,z2,z3) arrayfun(@(z1,z2,z3) func_s_array([z1;z2;z3]),z1,z2,z3);
end
end