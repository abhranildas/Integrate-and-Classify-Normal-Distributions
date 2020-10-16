function p_ray=prob_ray(init_sign_ray,z_ray,dim,side)

% function to compute probability slice along each ray
[f_big,f_small]=norm_rad_surv_split(z_ray,dim);
p_ray_big=init_sign_ray+1+init_sign_ray*sum((-1).^(1:length(z_ray)).*f_big);
p_ray_small=init_sign_ray*sum((-1).^(1:length(z_ray)).*f_small);
if ~exist('side','var')
    p_ray=p_ray_big+p_ray_small;
elseif strcmpi(side,'complement')
    p_ray=2-p_ray_big-p_ray_small;
end
end