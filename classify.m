function results=classify(dist_a,dist_b,varargin)

% parse inputs
p = inputParser;
validNum = @(x) isnumeric(x);
addRequired(p,'dist_a',validNum);
addRequired(p,'dist_b',validNum);
addParameter(p,'p_a',0.5, @(x) isnumeric(x) && isscalar(x) && (x > 0) && (x < 1));
addParameter(p,'custom_bd_coeffs',[]);
addParameter(p,'type','params', @(s) strcmp(s,'params') || strcmp(s,'obs'));
addParameter(p,'bPlot',1, @(x) (x == 0) || (x == 1));

parse(p,dist_a,dist_b,varargin{:});

if strcmp(p.Results.type,'params')
    mu_a=dist_a(:,1);
    vcov_a=dist_a(:,2:end);
    mu_b=dist_b(:,1);
    vcov_b=dist_b(:,2:end);
elseif strcmp(p.Results.type,'obs')
    mu_a=mean(dist_a)';
    vcov_a=cov(dist_a);
    mu_b=mean(dist_b)';
    vcov_b=cov(dist_b);    
end

% if input obs and p_a is not specified,
if strcmp(p.Results.type,'obs') && any(strcmp(p.UsingDefaults,'p_a'))
    % set p_a according to input category frequencies
    p_a=size(dist_a,1)/(size(dist_a,1)+size(dist_b,1));
else
    p_a=p.Results.p_a;
end    

p_b=1-p_a;
custom_bd_coeffs_a=p.Results.custom_bd_coeffs;
custom_bd_coeffs_b=custom_bd_coeffs_a;

bPlot=p.Results.bPlot;

dim=length(mu_a); % dimension

% if custom boundary is provided,
if ~isempty(custom_bd_coeffs_a)
    a2_c=custom_bd_coeffs_a.a2;
    a1_c=custom_bd_coeffs_a.a1;
    a0_c=custom_bd_coeffs_a.a0;
    % flip custom boundary sign for b
    custom_bd_coeffs_b.a2=-a2_c;
    custom_bd_coeffs_b.a1=-a1_c;
    custom_bd_coeffs_b.a0=-a0_c;
end

% approximate d' (considering equal vcov)
vcov_av=(vcov_a+vcov_b)/2; % assume each vcov = avg
d_gauss_aprx=sqrt((mu_a-mu_b)'/(vcov_av)*(mu_a-mu_b)); % d' with equal avg vcov & equal priors
%acc_opt_aprx=1-normcdf(d_aprx/2,'upper');

% Coefficients of optimal decision boundary:
a2_gauss=inv(vcov_b)-inv(vcov_a);
a1_gauss=2*(vcov_a\mu_a-vcov_b\mu_b);
a0_gauss=mu_b'/vcov_b*mu_b-mu_a'/vcov_a*mu_a+log((p_a/p_b)^2*det(vcov_b)/det(vcov_a));
bd_coeffs_gauss_opt=struct;
bd_coeffs_gauss_opt.a2=a2_gauss;
bd_coeffs_gauss_opt.a1=a1_gauss;
bd_coeffs_gauss_opt.a0=a0_gauss;

% compute optimal accuracy and d' for gaussians
[acc_gauss_opt_a]=accuracy_gauss(mu_a,vcov_a,mu_b,vcov_b,.5,[]);
[acc_gauss_opt_b]=accuracy_gauss(mu_b,vcov_b,mu_a,vcov_a,.5,[]);
acc_gauss_opt=mean([acc_gauss_opt_a,acc_gauss_opt_b]);
d_gauss=2*norminv(acc_gauss_opt);

% compute accuracy and boundary for each gaussian, and combine
[acc_gauss_a,dec_bd_a]=accuracy_gauss(mu_a,vcov_a,mu_b,vcov_b,p_a,custom_bd_coeffs_a);
[acc_gauss_b,dec_bd_b]=accuracy_gauss(mu_b,vcov_b,mu_a,vcov_a,p_b,custom_bd_coeffs_b);
acc_gauss=p_a*acc_gauss_a+p_b*acc_gauss_b;
bd_pts=[dec_bd_a,dec_bd_b];

% if no bdry found, dists are identical:
if ~numel(bd_pts)
    acc_gauss=0.5;
end

% if input is observations,
if strcmp(p.Results.type,'obs')
    % if custom boundary is provided,
    if ~isempty(custom_bd_coeffs_a)
        % compute observed accuracy with custom boundary
        [acc_obs,acc_obs_a,acc_obs_b]=accuracy_obs(dist_a,dist_b,custom_bd_coeffs_a);
    else
        %acc_obs_start=accuracy_obs(dist_a,dist_b,opt_bd_coeffs_gauss);
        % find boundary that optimizes observed accuracy
        fun = @(x)accuracy_obs_optimize(dim,x,dist_a,dist_b);
        x=fminsearch(fun,[a2_gauss(:); a1_gauss(:); a0_gauss]);
        if dim==1            
            a2_obs=x(1);
            a1_obs=x(2);
            a0_obs=x(3);
        elseif dim==2
            a2_obs=[x(1) x(2); x(3) x(4)];
            a1_obs=[x(5); x(6)];
            a0_obs=x(7);
        elseif dim==3
            a2_obs=[x(1) x(2) x(3); x(4) x(5) x(6); x(7) x(8) x(9)];
            a1_obs=[x(10); x(11); x(12)];
            a0_obs=x(13);
        end
        bd_coeffs_obs_opt=struct;
        bd_coeffs_obs_opt.a2=a2_obs;
        bd_coeffs_obs_opt.a1=a1_obs;
        bd_coeffs_obs_opt.a0=a0_obs;
        
        % optimal boundary points
        [~,bd_pts_obs_opt_a]=accuracy_gauss(mu_a,vcov_a,mu_b,vcov_b,p_a,bd_coeffs_obs_opt);
        [~,bd_pts_obs_opt_b]=accuracy_gauss(mu_b,vcov_b,mu_a,vcov_a,p_b,bd_coeffs_obs_opt);
        bd_pts_obs_opt=[bd_pts_obs_opt_a,bd_pts_obs_opt_b];
        
        % accuracy with optimal data boundary
        [acc_obs,acc_obs_a,acc_obs_b]=accuracy_obs(dist_a,dist_b,bd_coeffs_obs_opt);        
        if size(dist_a,1)==size(dist_b,1) % if both category frequencies same,
            d_obs=2*norminv(acc_obs); % obs accuracy can be used to compute obs d'
            results.d_obs=d_obs;
        end
    end
    results.acc_obs=acc_obs;
    results.acc_obs_a=acc_obs_a;
    results.acc_obs_b=acc_obs_b;
end

results.acc_gauss=acc_gauss;
results.acc_gauss_a=acc_gauss_a;
results.acc_gauss_b=acc_gauss_b;
results.d_gauss=d_gauss;
results.d_gauss_aprx=d_gauss_aprx;

if isempty(custom_bd_coeffs_a)    
    if strcmp(p.Results.type,'obs')
        results.bd_pts_obs_opt=bd_pts_obs_opt;
        results.bd_coeffs_obs_opt=bd_coeffs_obs_opt;        
    end
    results.bd_pts_gauss_opt=bd_pts;
else
    results.bd_pts_custom=bd_pts;
end
results.bd_coeffs_gauss_opt=bd_coeffs_gauss_opt;

% Plot:
if bPlot
    figure; hold on;
    % title    
    if isinf(d_gauss)
        title(sprintf("d'~%.2f",d_gauss_aprx),'color','red')
    else
        title(sprintf("d'=%.2f",d_gauss))
    end
    if strcmp(p.Results.type,'obs')
        if exist('d_obs','var')
            if ~isinf(d_obs)
                title(sprintf("d'=%.2f",d_obs))
            end
        end
    end
    
    %red=[196 52 38]/255; blue=[68 114 196]/255;
    
    if dim==1
        
        % plot data
        if strcmp(p.Results.type,'obs')
            [heights_a,edges_a]=histcounts(dist_a,'BinMethod','scott','normalization','pdf');
            histogram('BinCounts',heights_a*p_a,'BinEdges',edges_a,'facecolor','blue','facealpha',0.3,'edgecolor','blue','edgealpha',0.3);
            [heights_b,edges_b]=histcounts(dist_b,'BinMethod','scott','normalization','pdf');
            histogram('BinCounts',heights_b*p_b,'BinEdges',edges_b,'facecolor','red','facealpha',0.3,'edgecolor','red','edgealpha',0.3);
        end
        
        % plot gaussians
        x_a=linspace(mu_a-5*sqrt(vcov_a),mu_a+5*sqrt(vcov_a),100);
        y_a=p_a*normpdf(x_a,mu_a,sqrt(vcov_a));
        area(x_a,y_a,'facecolor','blue','facealpha',0.4,'edgecolor','blue','edgealpha',0.5,'linewidth',1)
        
        x_b=linspace(mu_b-5*sqrt(vcov_b),mu_b+5*sqrt(vcov_b),100);
        y_b=p_b*normpdf(x_b,mu_b,sqrt(vcov_b));
        area(x_b,y_b,'facecolor','red','facealpha',0.4,'edgecolor','red','edgealpha',0.5,'linewidth',1)
        
        % plot boundary
        for x=bd_pts
            line([x x],ylim,'color','k','linewidth',1)
        end
        if strcmp(p.Results.type,'obs')&& isempty(custom_bd_coeffs_a)
            for x=bd_pts_obs_opt
                line([x x],ylim,'color',.5*[1 1 1],'linewidth',1)
            end
        end
        
    elseif dim==2
        
        % plot data
        if strcmp(p.Results.type,'obs')
            scatter(dist_a(:,1),dist_a(:,2),4,'o','markerfacecolor','blue');
            scatter(dist_b(:,1),dist_b(:,2),4,'o','markerfacecolor','red');
        end
        
        % plot gaussians (error ellipses)
        th=0:.01:2*pi;
        z=[cos(th);sin(th)];        

        %plot(mu_a(1),mu_a(2),'.','markersize',20,'color','blue');
        C_a=chol(vcov_a,'lower');
        dist_a=C_a*z+repmat(mu_a,[1 length(th)]);
        
        %plot(mu_b(1),mu_b(2),'.','markersize',20,'color','red');
        C_b=chol(vcov_b,'lower');
        dist_b=C_b*z+repmat(mu_b,[1 length(th)]);
        
        plot(dist_a(1,:),dist_a(2,:),'-','color','blue')
        plot(dist_b(1,:),dist_b(2,:),'-','color','red')
        
        % plot boundary
        if numel(bd_pts)
            plot(bd_pts(1,:),bd_pts(2,:),'.k')
        end
        if strcmp(p.Results.type,'obs') && isempty(custom_bd_coeffs_a)
            plot(bd_pts_obs_opt(1,:),bd_pts_obs_opt(2,:),'.','color',.5*[1 1 1])
        end
        
        axis image
        %xlim([min(mu_a(1)-sig_a(1,1),mu_b(1)-sig_b(1,1)),max(mu_a(1)+sig_a(1,1),mu_b(1)+sig_b(1,1))])
        %ylim([min(mu_a(2)-sig_a(2,2),mu_b(2)-sig_b(2,2)),max(mu_a(2)+sig_a(2,2),mu_b(2)+sig_b(2,2))])
    
    elseif dim==3
        
        % plot data
        if strcmp(p.Results.type,'obs')
            plot3(dist_a(:,1),dist_a(:,2),dist_a(:,3),'.','color','blue','markersize',5);
            plot3(dist_b(:,1),dist_b(:,2),dist_b(:,3),'.','color','red','markersize',5);
        end
        
        % plot gaussians (error ellipsoids)
        [V,D] = eig(vcov_a); % V is the rotation matrix for the error ellipsoid.
        
        rot_angles=rad2deg(rotm2eul(V,'XYZ'));
        
        [x,y,z]=ellipsoid(mu_a(1),mu_a(2),mu_a(3),...
            sqrt(D(1,1)),sqrt(D(2,2)),sqrt(D(3,3)),20);
        
        ellipsoid_a=surf(x,y,z,'Facecolor','blue','Facealpha',.2,'Edgecolor','blue','Edgealpha',.3);
        
        rotate(ellipsoid_a,[1 0 0],rot_angles(1),mu_a)
        rotate(ellipsoid_a,[0 1 0],rot_angles(2),mu_a)
        rotate(ellipsoid_a,[0 0 1],rot_angles(3),mu_a)        

        [V,D] = eig(vcov_b); % V is the rotation matrix for the error ellipsoid.
        
        rot_angles=rad2deg(rotm2eul(V,'XYZ'));
        
        [x,y,z]=ellipsoid(mu_b(1),mu_b(2),mu_b(3),...
            sqrt(D(1,1)),sqrt(D(2,2)),sqrt(D(3,3)),20);
        
        ellipsoid_b=surf(x,y,z,'Facecolor','red','Facealpha',.2,'Edgecolor','red','Edgealpha',.3);
        
        rotate(ellipsoid_b,[1 0 0],rot_angles(1),mu_b)
        rotate(ellipsoid_b,[0 1 0],rot_angles(2),mu_b)
        rotate(ellipsoid_b,[0 0 1],rot_angles(3),mu_b)
        
        % plot boundary       
        if numel(bd_pts)
            %plot3(bd_pts(1,:),bd_pts(2,:),bd_pts(3,:),'.k','markersize',5);
            [t]=MyCrustOpen(bd_pts');
            trisurf(t,bd_pts(1,:),bd_pts(2,:),bd_pts(3,:),'facecolor',.5*[1 1 1],'facealpha',.1,'edgecolor','k','edgealpha',.2)
        end
        if strcmp(p.Results.type,'obs')
            plot3(bd_pts_obs_opt(1,:),bd_pts_obs_opt(2,:),bd_pts_obs_opt(3,:),'.','color',.5*[1 1 1],'markersize',5)
        end
        
        axis image
        grid on
        xlim([min(mu_a(1)-vcov_a(1,1),mu_b(1)-vcov_b(1,1)),max(mu_a(1)+vcov_a(1,1),mu_b(1)+vcov_b(1,1))])
        ylim([min(mu_a(2)-vcov_a(2,2),mu_b(2)-vcov_b(2,2)),max(mu_a(2)+vcov_a(2,2),mu_b(2)+vcov_b(2,2))])
        zlim([min(mu_a(3)-vcov_a(3,3),mu_b(3)-vcov_b(3,3)),max(mu_a(3)+vcov_a(3,3),mu_b(3)+vcov_b(3,3))])

        xlabel('X')
        ylabel('Y')
        zlabel('Z')       

    end
    hold off
end