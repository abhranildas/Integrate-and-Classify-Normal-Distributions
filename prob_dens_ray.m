function f_ray=prob_dens_ray(n_ray,z_ray,standard_gradf)

    % computes the probability density of a function of a multinormal along
    % a ray

    % parse inputs
    parser=inputParser;
    parser.KeepUnmatched=true;
    addRequired(parser,'n_ray',@isnumeric); % ray vector
    addRequired(parser,'z_ray',@isnumeric); % function roots on the ray
    addRequired(parser,'standard_gradf',@(x) isa(x,'function_handle')); % gradient of standardized function

    parse(parser,n_ray,z_ray,standard_gradf);
    dim=size(n_ray,1);
    
    % slope of function at roots on the ray
    fprime=arrayfun(@(z) dot(standard_gradf(n_ray*z),n_ray),z_ray);

    f_ray=sum(phi_ray(z_ray,dim)./abs(fprime));

end