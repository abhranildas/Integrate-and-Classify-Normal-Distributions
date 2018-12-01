function [dPrime,err_rate,dec_bd,dec_bd_sign,bPlot]=classification_error(mu_a,sig_a,mu_b,sig_b,p_a,bPlot)

d=length(mu_a);

[err_rate_a,dec_bd_a,dec_bd_a_sign]=mass_outside(mu_a,sig_a,mu_b,sig_b,p_a);
[err_rate_b,dec_bd_b,dec_bd_b_sign]=mass_outside(mu_b,sig_b,mu_a,sig_a,1-p_a);

dec_bd=[dec_bd_a,dec_bd_b];
dec_bd_sign=[dec_bd_a_sign,dec_bd_b_sign];

err_rate=(err_rate_a+err_rate_b)/2;
dPrime=-2*norminv(err_rate);

% Plot distributions and decision boundary:
if bPlot
    figure; hold on
    red=[1 0 0]; green=[0 1 0]; blue=[0 0 1];
    
    if d==2
        th=0:.01:2*pi;
        z=[cos(th);sin(th)];
        
        % Error ellipse a
        C_a=chol(sig_a,'lower');
        dist_a=C_a*z+repmat(mu_a,[1 length(th)]);
        
        % Error ellipse b
        C_b=chol(sig_b,'lower');
        dist_b=C_b*z+repmat(mu_b,[1 length(th)]);
        
        
        plot(dist_a(1,:),dist_a(2,:),'-b')
        plot(dist_b(1,:),dist_b(2,:),'-r')
        
        plot(dec_bd(1,:),dec_bd(2,:),'.k')
        
        axis image
        xlim([min(mu_a(1)-sig_a(1,1),mu_b(1)-sig_b(1,1)),max(mu_a(1)+sig_a(1,1),mu_b(1)+sig_b(1,1))])
        ylim([min(mu_a(2)-sig_a(2,2),mu_b(2)-sig_b(2,2)),max(mu_a(2)+sig_a(2,2),mu_b(2)+sig_b(2,2))])
    
    elseif d==3
        
        % Error ellipsoid a
        [V,D] = eig(sig_a); % V is the rotation matrix for the error ellipsoid.
        
        rot_angles=rad2deg(rotm2eul(V,'XYZ'));
        
        [x,y,z]=ellipsoid(mu_a(1),mu_a(2),mu_a(3),...
            sqrt(D(1,1)),sqrt(D(2,2)),sqrt(D(3,3)),20);
        
        ellipsoid_a=surf(x,y,z,'Facecolor',blue,'Facealpha',.4,'Edgecolor',blue,'Edgealpha',.6);
        
        rotate(ellipsoid_a,[1 0 0],rot_angles(1),mu_a)
        rotate(ellipsoid_a,[0 1 0],rot_angles(2),mu_a)
        rotate(ellipsoid_a,[0 0 1],rot_angles(3),mu_a)
        
        % Error ellipsoid b
        [V,D] = eig(sig_b); % V is the rotation matrix for the error ellipsoid.
        
        rot_angles=rad2deg(rotm2eul(V,'XYZ'));
        
        [x,y,z]=ellipsoid(mu_b(1),mu_b(2),mu_b(3),...
            sqrt(D(1,1)),sqrt(D(2,2)),sqrt(D(3,3)),20);
        
        ellipsoid_b=surf(x,y,z,'Facecolor',red,'Facealpha',.4,'Edgecolor',red,'Edgealpha',.6);
        
        rotate(ellipsoid_b,[1 0 0],rot_angles(1),mu_b)
        rotate(ellipsoid_b,[0 1 0],rot_angles(2),mu_b)
        rotate(ellipsoid_b,[0 0 1],rot_angles(3),mu_b)
        
        % Decision surface        
        scatter3(dec_bd(1,:),dec_bd(2,:),dec_bd(3,:),15,'.k');
        %[t]=MyCrustOpen(dec_bd_a');
        %trisurf(t,dec_bd_b(1,:),dec_bd_b(2,:),dec_bd_b(3,:),'facecolor',.5*[1 1 1],'facealpha',.1,'edgecolor','k','edgealpha',.2)
        
        
        axis equal
        grid on
        xlim([min(mu_a(1)-sig_a(1,1),mu_b(1)-sig_b(1,1)),max(mu_a(1)+sig_a(1,1),mu_b(1)+sig_b(1,1))])
        ylim([min(mu_a(2)-sig_a(2,2),mu_b(2)-sig_b(2,2)),max(mu_a(2)+sig_a(2,2),mu_b(2)+sig_b(2,2))])
        zlim([min(mu_a(3)-sig_a(3,3),mu_b(3)-sig_b(3,3)),max(mu_a(3)+sig_a(3,3),mu_b(3)+sig_b(3,3))])

        xlabel('X')
        ylabel('Y')
        zlabel('Z')       

    end
    hold off
end