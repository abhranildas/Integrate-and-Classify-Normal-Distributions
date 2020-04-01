function [f_b,f_s]=norm_rad_surv_split(z,dim)
z_c=chi2inv(0.5,dim);
f_b=((z<=-z_c)+(z<z_c))/2;

    function f_s_one=f_small(z_one)
        if abs(z_one)<z_c
            f_s_one=chi2cdf(z_one^2,dim)/2;
        else
            f_s_one=chi2cdf(z_one^2,dim,'upper')/2;
        end
        if (z_one<=-z_c)||((z_one>0)&&(z_one<z_c))
            f_s_one=-f_s_one;
        end
    end
f_s=arrayfun(@f_small,z);
end