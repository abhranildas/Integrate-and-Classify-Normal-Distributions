function plot_normal(mu,v,bd_pts,prior,color)
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
    colors=colororder;
    color=colors(1,:);
end

dim=length(mu);

if dim==1
    % plot normal
    x=linspace(mu-5*sqrt(v),mu+5*sqrt(v),100);
    y=prior*normpdf(x,mu,sqrt(v));
    area(x,y,'facecolor',color,'facealpha',0.4,'edgecolor',color,'edgealpha',0.5,'linewidth',1)
    
    % plot boundary
    hold on
    for x=bd_pts
        xline(x,'linewidth',1);
    end
    
elseif dim==2
    % plot normal (error ellipse)
    th=0:.01:2*pi;
    z=[cos(th);sin(th)];
    C=chol(v,'lower');
    dist=C*z+repmat(mu,[1 length(th)]);
    fill(dist(1,:),dist(2,:),color,'linestyle','none','facealpha',.5);
    
    % plot boundary
    hold on
    if numel(bd_pts)
        plot(bd_pts(1,:),bd_pts(2,:),'.k','markersize',3)
    end
    
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
    hold on
    if numel(bd_pts)
        plot3(bd_pts(1,:),bd_pts(2,:),bd_pts(3,:),'.k','markersize',1);
    end
    grid on
end
hold off