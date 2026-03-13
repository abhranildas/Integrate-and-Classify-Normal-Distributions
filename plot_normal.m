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

    % [R,D] = eig(C); % R is the rotation matrix for the error ellipsoid.
    % rot_angles=rad2deg(rotm2eul(R,'XYZ'));
    % [x,y,z]=ellipsoid(mu(1),mu(2),mu(3),D(1,1),D(2,2),D(3,3),20);
    % err_ellipsoid=surf(x,y,z,'Facecolor',plot_color,'Facealpha',.2,'Edgecolor',plot_color,'Edgealpha',.3);
    % rotate(err_ellipsoid,[1 0 0],rot_angles(1),mu)
    % rotate(err_ellipsoid,[0 1 0],rot_angles(2),mu)
    % rotate(err_ellipsoid,[0 0 1],rot_angles(3),mu)


    % v: 3x3 covariance matrix
    % mu: 3x1 or 1x3 mean (center)
    % plot_color: your desired color, e.g. 'r' or [0 0.447 0.741]


    % Eigen-decomposition of covariance
    %    Columns of R are principal directions
    [R, Lambda] = eig(v);      % v * R = R * Lambda

    % Radii = sqrt of eigenvalues (std devs along principal axes)
    radii = sqrt(diag(Lambda));

    % Ellipsoid in its own principal-axis coordinate system
    %    Centered at origin for now
    [xe, ye, ze] = ellipsoid(0, 0, 0, radii(1), radii(2), radii(3), 20);

    % Flatten, rotate, then translate
    pts = [xe(:)'; ye(:)'; ze(:)'];  % 3 x N

    % Apply rotation given by eigenvectors (columns of R)
    pts_rot = R * pts;              % 3 x N

    % Add the mean to move to correct center
    pts_rot(1,:) = pts_rot(1,:) + mu(1);
    pts_rot(2,:) = pts_rot(2,:) + mu(2);
    pts_rot(3,:) = pts_rot(3,:) + mu(3);

    % 6) Reshape back to grid form
    x = reshape(pts_rot(1,:), size(xe));
    y = reshape(pts_rot(2,:), size(xe));
    z = reshape(pts_rot(3,:), size(xe));

    % Plot
    surf(x, y, z,'FaceColor',plot_color,'FaceAlpha', 0.2,'EdgeColor',plot_color,'EdgeAlpha', 0.3);
    view([0 45])
end

% set axis limits
if dim>1
    ax_curr=reshape(axis,2,[])';
    ax_this=mu+3*sqrt(diag(v)).*[-1 1];
    ax_combined=[min([ax_curr(:,1) ax_this(:,1)],[],2),max([ax_curr(:,2) ax_this(:,2)],[],2)];
    axis(reshape(ax_combined',1,[]));
end

if ~holdon
    hold off
end
