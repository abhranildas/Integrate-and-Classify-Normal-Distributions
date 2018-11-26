function [d,err_rate,dec_bd_raw,dec_bd_sign,bPlot]=d_prime(mu_a,sig_a,mu_b,sig_b,bPlot)
    
[err_rate_a,dec_bd_a]=mass_outside(mu_a,sig_a,mu_b,sig_b);
[err_rate_b,dec_bd_b]=mass_outside(mu_b,sig_b,mu_a,sig_a);

err_rate=(err_rate_a+err_rate_b)/2;
d=-2*norminv(err_rate);

% plot distributions and decision bound:
if bPlot
    th=0:.01:2*pi;
    z=[cos(th);sin(th)];
    
    C_a=chol(sig_a,'lower');
    dist_a=C_a*z+repmat(mu_a,[1 length(th)]);
    
    C_b=chol(sig_b,'lower');
    dist_b=C_b*z+repmat(mu_b,[1 length(th)]);
    
    figure; hold on
    plot(dist_a(1,:),dist_a(2,:),'-b')
    plot(dist_b(1,:),dist_b(2,:),'-r')
    
    plot(dec_bd_a(1,:),dec_bd_a(2,:),'.k')
    plot(dec_bd_b(1,:),dec_bd_b(2,:),'.k')
    
    axis image
    xlim([min(mu_a(1)-sig_a(1,1),mu_b(1)-sig_b(1,1)),max(mu_a(1)+sig_a(1,1),mu_b(1)+sig_b(1,1))])
    ylim([min(mu_a(2)-sig_a(2,2),mu_b(2)-sig_b(2,2)),max(mu_a(2)+sig_a(2,2),mu_b(2)+sig_b(2,2))])
end