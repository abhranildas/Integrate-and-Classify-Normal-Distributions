function p_ray=prob_ray_sym(init_sign_ray,z_ray,dim,side)

    % function to compute probability slice along each ray symbolically, for
    % tiny probabilities

    sgns=sym(init_sign_ray);
    ds=sym(dim);

    % complementary cdf on ray using symbol-friendly igamma
    % first use igamma to construct regularized gamma, using which
    % construct chi distribution cdf, using which construct phibar.

    function p=Phibar(z,ds)
        % chi CDF is regularized gamma, built from symbol-friendly igamma
        chiCDF=1-igamma(ds/2,z^2/2)/gamma(ds/2);

        % Phibar from chicdf
        p=(1-sign(z)*chiCDF)/2;
    end

    Phibars=arrayfun(@(z) Phibar(z,ds),z_ray);

    p_ray=sgns+1+2*sgns*sum((-1).^(1:length(z_ray)).*Phibars);

    if exist('side','var') && strcmpi(side,'lower')
        p_ray=2-p_ray;
    end

end