function [f_big,f_small]=norm_rad_surv_split(z,dim)
% returns the normal radial survival function, split into a big chunk, which is 
% 0, +/-2 or +/-1, and a small part, which is a tail CDF. This is to
% prevent the small part vanishing when the two are added together due to machine imprecision.

z_c=chi2inv(0.5,dim); % above this criterion, we'll use chi2cdf('upper').
f_big=(z<=-z_c)+(z<z_c); % 2 when z<=-z_c, 1 when -z_c<z<z_c, 0 when z>=z_c

f_small=chi2cdf(z.^2,dim);
f_small_upper=chi2cdf(z.^2,dim,'upper');
f_small(abs(z)>=z_c)=f_small_upper(abs(z)>=z_c);
f_small=f_small.*(1-2*((z<=-z_c)|((z>0)&(z<z_c)))); % if z<=-z_c or 0<z<z_c, invert sign
end