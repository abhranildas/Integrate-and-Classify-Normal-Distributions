function plot_normal(mu,v,prior,plot_color)
    % Plot a normal distribution and a boundary.
    %
    % Credits:
    %   Abhranil Das <abhranil.das@utexas.edu>
    %	Wilson S Geisler
    %	Center for Perceptual Systems, University of Texas at Austin
    % If you use this code, please cite:
    %   A new method to compute classification error
    %   https://jov.arvojournals.org/article.aspx?articleid=2750251

    if ~exist('prior','var')
        prior=1;
    end

    if ~exist('plot_color','var')
        colors=colororder;
        plot_color=colors(1,:);
    end

    dim=length(mu);

    C=sqrtm(v);

    holdon=ishold;

    if dim==1
        % plot normal
        fplot(@(x) prior*normpdf(x,mu,C), 'color',plot_color,'linewidth',1);
        hold on
        x=linspace(mu-5*C,mu+5*C,100);
        y=prior*normpdf(x,mu,C);
        area(x,y,'facecolor',plot_color,'facealpha',0.4,'edgecolor','none')

    elseif dim==2
        % plot normal (error ellipse)
        th=linspace(0,2*pi,100);
        z=[cos(th);sin(th)];
        dist=C*z+mu;
        fill(dist(1,:),dist(2,:),plot_color,'edgecolor','none','facealpha',.5);

    elseif dim==3
        % plot normal (error ellipsoid)
        grid on
        [R,D] = eig(C); % R is the rotation matrix for the error ellipsoid.
        rot_angles=rad2deg(rotm2eul(R,'XYZ'));
        [x,y,z]=ellipsoid(mu(1),mu(2),mu(3),D(1,1),D(2,2),D(3,3),20);
        err_ellipsoid=surf(x,y,z,'Facecolor',plot_color,'Facealpha',.2,'Edgecolor',plot_color,'Edgealpha',.3);
        rotate(err_ellipsoid,[1 0 0],rot_angles(1),mu)
        rotate(err_ellipsoid,[0 1 0],rot_angles(2),mu)
        rotate(err_ellipsoid,[0 0 1],rot_angles(3),mu)
        view([0 45])
    end

    % set axis limits
    if dim>1
        ax_curr=reshape(axis,2,[])';
        ax_this=mu+3*sqrt(diag(v)).*[-1 1];
        ax_combined=[min([ax_curr(:,1) ax_this(:,1)],[],2),max([ax_curr(:,2) ax_this(:,2)],[],2)];
        axis image; axis(reshape(ax_combined',1,[]));
    end

    if ~holdon
        hold off
    end
