function plot_normal(mu,v,bd_pts,p_prior,color)
% Plot a normal distribution and a boundary.
%
% Credits:
%   Abhranil Das <abhranil.das@utexas.edu>
%	Wilson S Geisler
%	Center for Perceptual Systems, University of Texas at Austin
% If you use this code, please cite:
%   A new method to compute classification error
%   https://jov.arvojournals.org/article.aspx?articleid=2750251

if ~exist('color','var')
    color='blue';
end

hold on

dim=length(mu);

if dim==1
    % plot normal
    x=linspace(mu-5*sqrt(v),mu+5*sqrt(v),100);
    y=p_prior*normpdf(x,mu,sqrt(v));
    area(x,y,'facecolor',color,'facealpha',0.4,'edgecolor',color,'edgealpha',0.5,'linewidth',1)
    
    % plot boundary
    for x=bd_pts
        line([x x],ylim,'color','k','linewidth',1)
    end
    
elseif dim==2
    % plot normal (error ellipse)
    th=0:.01:2*pi;
    z=[cos(th);sin(th)];
    %plot(mu_a(1),mu_a(2),'.','markersize',20,'color','blue');
    C=chol(v,'lower');
    dist=C*z+repmat(mu,[1 length(th)]);
    fill(dist(1,:),dist(2,:),color,'facealpha',.5);
    
    % plot boundary
    if numel(bd_pts)
        plot(bd_pts(1,:),bd_pts(2,:),'.k')
    end
    
    %axis image
    %xlim([min(mu_a(1)-sig_a(1,1),mu_b(1)-sig_b(1,1)),max(mu_a(1)+sig_a(1,1),mu_b(1)+sig_b(1,1))])
    %ylim([min(mu_a(2)-sig_a(2,2),mu_b(2)-sig_b(2,2)),max(mu_a(2)+sig_a(2,2),mu_b(2)+sig_b(2,2))])
    
elseif dim==3    
    % plot normal (error ellipsoid)
    [V,D] = eig(v); % V is the rotation matrix for the error ellipsoid.    
    rot_angles=rad2deg(rotm2eul(V,'XYZ'));    
    [x,y,z]=ellipsoid(mu(1),mu(2),mu(3),sqrt(D(1,1)),sqrt(D(2,2)),sqrt(D(3,3)),20);
    err_ellipsoid=surf(x,y,z,'Facecolor',color,'Facealpha',.2,'Edgecolor',color,'Edgealpha',.3);
    rotate(err_ellipsoid,[1 0 0],rot_angles(1),mu)
    rotate(err_ellipsoid,[0 1 0],rot_angles(2),mu)
    rotate(err_ellipsoid,[0 0 1],rot_angles(3),mu)
    
    % plot boundary
    if numel(bd_pts)
        plot3(bd_pts(1,:),bd_pts(2,:),bd_pts(3,:),'.k','markersize',5);
        %[t]=MyCrustOpen(bd_pts');
        %trisurf(t,bd_pts(1,:),bd_pts(2,:),bd_pts(3,:),'facecolor',.5*[1 1 1],'facealpha',.1,'edgecolor','k','edgealpha',.2)
    end
    %axis image
    grid on
    xlabel x
    ylabel y
    zlabel z    
end
hold off