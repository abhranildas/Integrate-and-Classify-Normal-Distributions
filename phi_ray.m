function f=phi_ray(z,dim)
    % standard multinormal density function along a ray
    f=abs(z).*chi2pdf(z.^2,dim);
end