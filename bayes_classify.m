function results=bayes_classify(dist_a,dist_b,varargin)
% Classification between two distributions.
% Credits:
%   Abhranil Das <abhranil.das@utexas.edu>
%	R Calen Walshe
%	Wilson S Geisler
%	Center for Perceptual Systems, University of Texas at Austin
% If you use this code, please cite:
%   A new method to compute classification error
%   https://jov.arvojournals.org/article.aspx?articleid=2750251

% parse inputs
p = inputParser;
addRequired(p,'dist_a',@(x) isnumeric(x));
addRequired(p,'dist_b',@(x) isnumeric(x));
addParameter(p,'p_a',0.5, @(x) isnumeric(x) && isscalar(x) && (x > 0) && (x < 1));
addParameter(p,'vals',eye(2), @(x) isnumeric(x) && ismatrix(x));
addParameter(p,'custom_bd_coeffs',[]);
addParameter(p,'type','params', @(s) strcmp(s,'params') || strcmp(s,'obs'));
addParameter(p,'bPlot',1, @(x) (x == 0) || (x == 1));

parse(p,dist_a,dist_b,varargin{:});

if strcmp(p.Results.type,'params')
    mu_a=dist_a(:,1);
    v_a=dist_a(:,2:end);
    mu_b=dist_b(:,1);
    v_b=dist_b(:,2:end);
elseif strcmp(p.Results.type,'obs')
    mu_a=mean(dist_a)';
    v_a=cov(dist_a);
    mu_b=mean(dist_b)';
    v_b=cov(dist_b);    
end

% if input obs and p_a is not specified,
if strcmp(p.Results.type,'obs') && any(strcmp(p.UsingDefaults,'p_a'))
    % set p_a according to input sample sizes
    p_a=size(dist_a,1)/(size(dist_a,1)+size(dist_b,1));
else
    p_a=p.Results.p_a;
end    

p_b=1-p_a;

vals=p.Results.vals;

custom_bd_coeffs_a=p.Results.custom_bd_coeffs;
custom_bd_coeffs_b=custom_bd_coeffs_a;

bPlot=p.Results.bPlot;

dim=length(mu_a); % dimension

% if custom boundary is provided,
if ~isempty(custom_bd_coeffs_a)
    % flip custom boundary sign for b
    custom_bd_coeffs_b.a2=-custom_bd_coeffs_a.a2;
    custom_bd_coeffs_b.a1=-custom_bd_coeffs_a.a1;
    custom_bd_coeffs_b.a0=-custom_bd_coeffs_a.a0;
end

% approximate d' (considering equal vcov)
vcov_av=(v_a+v_b)/2; % assume each vcov = avg
d_gauss_aprx=sqrt((mu_a-mu_b)'/(vcov_av)*(mu_a-mu_b)); % d' with equal avg vcov & equal priors
results.d_gauss_aprx=d_gauss_aprx;
% acc_opt_gauss_aprx=1-normcdf(d_gauss_aprx/2,'upper');
% results.acc_opt_gauss_aprx=acc_opt_gauss_aprx;

% Chernoff bound for error and d'
fun = @(b)chernoff_bound(b,mu_a,v_a,mu_b,v_b,p_a);
[~,k]=fminbnd(fun,0,1);
log_err_gauss_max=k;
results.log_err_gauss_max=log_err_gauss_max;

fun = @(b)chernoff_bound(b,mu_a,v_a,mu_b,v_b,0.5);
[~,k]=fminbnd(fun,0,1);
acc_opt_min_chernoff=1-10^k;
d_gauss_min=2*norminv(acc_opt_min_chernoff);
results.d_gauss_min=d_gauss_min;

% Coefficients of optimal decision boundary:
a2_gauss=inv(v_b)-inv(v_a);
a1_gauss=2*(v_a\mu_a-v_b\mu_b);
a0_gauss=mu_b'/v_b*mu_b-mu_a'/v_a*mu_a+log((((vals(1,1)-vals(1,2))*p_a)/((vals(2,2)-vals(2,1))*p_b))^2*det(v_b)/det(v_a));
bd_coeffs_gauss_opt_a=struct;
bd_coeffs_gauss_opt_a.a2=a2_gauss;
bd_coeffs_gauss_opt_a.a1=a1_gauss;
bd_coeffs_gauss_opt_a.a0=a0_gauss;
results.bd_coeffs_gauss_opt=bd_coeffs_gauss_opt_a;

if dim<=3
    % compute optimal accuracy and d' for gaussians
    [acc_gauss_opt_a]=accuracy_gauss(mu_a,v_a,mu_b,v_b,.5,eye(2),[]);
    [acc_gauss_opt_b]=accuracy_gauss(mu_b,v_b,mu_a,v_a,.5,eye(2),[]);
    acc_gauss_opt=mean([acc_gauss_opt_a,acc_gauss_opt_b]);
    d_gauss=2*norminv(acc_gauss_opt);
    
    % compute accuracy and boundary for each gaussian, and combine
    [acc_gauss_a,bd_pts_gauss_opt_a]=accuracy_gauss(mu_a,v_a,mu_b,v_b,p_a,vals,custom_bd_coeffs_a);
    [acc_gauss_b,bd_pts_gauss_opt_b]=accuracy_gauss(mu_b,v_b,mu_a,v_a,p_b,rot90(vals,2),custom_bd_coeffs_b);
    acc_gauss=p_a*acc_gauss_a+p_b*acc_gauss_b;
    bd_pts_gauss_opt=[bd_pts_gauss_opt_a,bd_pts_gauss_opt_b];
    
    if isempty(custom_bd_coeffs_a)
        results.bd_pts_gauss_opt=bd_pts_gauss_opt;
    else
        results.bd_pts_custom=bd_pts_gauss_opt;
    end
    
    % if no bdry found, dists are identical:
    if ~numel(bd_pts_gauss_opt)
        acc_gauss=0.5;
    end
    
    results.errs_gauss=[1-acc_gauss, 1-acc_gauss_a, 1-acc_gauss_b];
%     results.acc_gauss_a=acc_gauss_a;
%     results.acc_gauss_b=acc_gauss_b;
    results.d_gauss=d_gauss;  
    
    if ~isequal(vals,eye(2)) % if outcome values are supplied
        results.outcome_vals_gauss=[vals(1,1)*p_a*acc_gauss_a, vals(1,2)*p_a*(1-acc_gauss_a);...
            vals(2,1)*p_b*(1-acc_gauss_b), vals(2,2)*p_b*acc_gauss_b];
        results.val_ex_gauss=(vals(1,1)-vals(1,2))*p_a*acc_gauss_a + ...
            (vals(2,2)-vals(2,1))*p_b*acc_gauss_b + vals(1,2)*p_a + vals(2,1)*p_b;
    end
end

% if input is observations,
if strcmp(p.Results.type,'obs')
    % if custom boundary is provided,
    if ~isempty(custom_bd_coeffs_a)
        % compute accuracy and expected value with custom boundary
        [acc_obs,acc_obs_a,acc_obs_b,outcome_counts_obs]=val_obs(dist_a,dist_b,custom_bd_coeffs_a,eye(size(dist_a,2)));
        [val_ex_obs, val_ex_obs_a, val_ex_obs_b, outcome_vals_obs]=val_obs(dist_a,dist_b,custom_bd_coeffs_a,vals);
    else
        %acc_obs_start=accuracy_obs(dist_a,dist_b,opt_bd_coeffs_gauss);
        % find boundary that optimizes expected value / accuracy
        fun = @(x)-val_obs_flat(dim,x,dist_a,dist_b,vals);
        x=fminsearch(fun,[a2_gauss(:); a1_gauss(:); a0_gauss],optimset('Display','iter'));

        a2_obs=reshape(x(1:dim^2),[dim dim])';
        a1_obs=x(dim^2+1:dim^2+dim);
        a0_obs=x(end);

        bd_coeffs_obs_opt_a=struct;
        bd_coeffs_obs_opt_a.a2=a2_obs;
        bd_coeffs_obs_opt_a.a1=a1_obs;
        bd_coeffs_obs_opt_a.a0=a0_obs;
        results.bd_coeffs_obs_opt=bd_coeffs_obs_opt_a;
        
        % flip boundary sign for b
        bd_coeffs_obs_opt_b.a2=-bd_coeffs_obs_opt_a.a2;
        bd_coeffs_obs_opt_b.a1=-bd_coeffs_obs_opt_a.a1;
        bd_coeffs_obs_opt_b.a0=-bd_coeffs_obs_opt_a.a0;
        
        if dim<=3
            % optimal boundary points
            [~,bd_pts_obs_opt_a]=accuracy_gauss(mu_a,v_a,[],[],[],[],bd_coeffs_obs_opt_a);
            [~,bd_pts_obs_opt_b]=accuracy_gauss(mu_b,v_b,[],[],[],[],bd_coeffs_obs_opt_b);
            bd_pts_obs_opt=[bd_pts_obs_opt_a,bd_pts_obs_opt_b];
            results.bd_pts_obs_opt=bd_pts_obs_opt;
        end
        
        % accuracy with optimal data boundary
        [acc_obs,acc_obs_a,acc_obs_b,outcome_counts_obs]=val_obs(dist_a,dist_b,bd_coeffs_obs_opt_a,eye(2));
        
        % expected value with optimal data boundary
        [val_ex_obs,~,~,outcome_vals_obs]=val_obs(dist_a,dist_b,bd_coeffs_obs_opt_a,vals);
        
        if size(dist_a,1)==size(dist_b,1) % if both category frequencies same,
            d_obs=2*norminv(acc_obs); % obs accuracy can be used to compute obs d'
            results.d_obs=d_obs;
        end
    end
    results.errs_obs=[1-acc_obs,1-acc_obs_a,1-acc_obs_b];
    results.outcome_counts_obs=outcome_counts_obs;
    if ~isequal(vals,eye(2)) % if outcome values are supplied
        results.val_ex_obs=val_ex_obs;
        results.outcome_vals_obs=outcome_vals_obs;
    end
end


% Plot:
if bPlot && dim<=3
    figure; hold on;
    % title
    if strcmp(p.Results.type,'params')
        title(sprintf("error = %g",1-acc_gauss(1)))
    else
        title(sprintf("error = %g",1-acc_obs(1)))
    end
    
    if dim==1
        
        % plot data
        if strcmp(p.Results.type,'obs')
            [heights_a,edges_a]=histcounts(dist_a,'BinMethod','scott','normalization','pdf');
            histogram('BinCounts',heights_a*p_a,'BinEdges',edges_a,'facecolor','blue','facealpha',0.3,'edgecolor','blue','edgealpha',0.3);
            [heights_b,edges_b]=histcounts(dist_b,'BinMethod','scott','normalization','pdf');
            histogram('BinCounts',heights_b*p_b,'BinEdges',edges_b,'facecolor','red','facealpha',0.3,'edgecolor','red','edgealpha',0.3);
        end
        
        % plot gaussians
        x_a=linspace(mu_a-5*sqrt(v_a),mu_a+5*sqrt(v_a),100);
        y_a=p_a*normpdf(x_a,mu_a,sqrt(v_a));
        area(x_a,y_a,'facecolor','blue','facealpha',0.4,'edgecolor','blue','edgealpha',0.5,'linewidth',1)
        
        x_b=linspace(mu_b-5*sqrt(v_b),mu_b+5*sqrt(v_b),100);
        y_b=p_b*normpdf(x_b,mu_b,sqrt(v_b));
        area(x_b,y_b,'facecolor','red','facealpha',0.4,'edgecolor','red','edgealpha',0.5,'linewidth',1)
        
        % plot boundary
        for x=bd_pts_gauss_opt
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
        C_a=chol(v_a,'lower');
        dist_a=C_a*z+repmat(mu_a,[1 length(th)]);
        
        %plot(mu_b(1),mu_b(2),'.','markersize',20,'color','red');
        C_b=chol(v_b,'lower');
        dist_b=C_b*z+repmat(mu_b,[1 length(th)]);
        
        plot(dist_a(1,:),dist_a(2,:),'-','color','blue')
        plot(dist_b(1,:),dist_b(2,:),'-','color','red')
        
        % plot boundary
        if numel(bd_pts_gauss_opt)
            plot(bd_pts_gauss_opt(1,:),bd_pts_gauss_opt(2,:),'.k')
        end
        if strcmp(p.Results.type,'obs') && isempty(custom_bd_coeffs_a)
            plot(bd_pts_obs_opt(1,:),bd_pts_obs_opt(2,:),'.','color',.5*[1 1 1])
        end
        
        %axis image
        %xlim([min(mu_a(1)-sig_a(1,1),mu_b(1)-sig_b(1,1)),max(mu_a(1)+sig_a(1,1),mu_b(1)+sig_b(1,1))])
        %ylim([min(mu_a(2)-sig_a(2,2),mu_b(2)-sig_b(2,2)),max(mu_a(2)+sig_a(2,2),mu_b(2)+sig_b(2,2))])
    
    elseif dim==3
        
        % plot data
        if strcmp(p.Results.type,'obs')
            plot3(dist_a(:,1),dist_a(:,2),dist_a(:,3),'.','color','blue','markersize',5);
            plot3(dist_b(:,1),dist_b(:,2),dist_b(:,3),'.','color','red','markersize',5);
        end
        
        % plot gaussians (error ellipsoids)
        [V,D] = eig(v_a); % V is the rotation matrix for the error ellipsoid.
        
        rot_angles=rad2deg(rotm2eul(V,'XYZ'));
        
        [x,y,z]=ellipsoid(mu_a(1),mu_a(2),mu_a(3),...
            sqrt(D(1,1)),sqrt(D(2,2)),sqrt(D(3,3)),20);
        
        ellipsoid_a=surf(x,y,z,'Facecolor','blue','Facealpha',.2,'Edgecolor','blue','Edgealpha',.3);
        
        rotate(ellipsoid_a,[1 0 0],rot_angles(1),mu_a)
        rotate(ellipsoid_a,[0 1 0],rot_angles(2),mu_a)
        rotate(ellipsoid_a,[0 0 1],rot_angles(3),mu_a)        

        [V,D] = eig(v_b); % V is the rotation matrix for the error ellipsoid.
        
        rot_angles=rad2deg(rotm2eul(V,'XYZ'));
        
        [x,y,z]=ellipsoid(mu_b(1),mu_b(2),mu_b(3),...
            sqrt(D(1,1)),sqrt(D(2,2)),sqrt(D(3,3)),20);
        
        ellipsoid_b=surf(x,y,z,'Facecolor','red','Facealpha',.2,'Edgecolor','red','Edgealpha',.3);
        
        rotate(ellipsoid_b,[1 0 0],rot_angles(1),mu_b)
        rotate(ellipsoid_b,[0 1 0],rot_angles(2),mu_b)
        rotate(ellipsoid_b,[0 0 1],rot_angles(3),mu_b)
        
        % plot boundary       
        if numel(bd_pts_gauss_opt)
            plot3(bd_pts_gauss_opt(1,:),bd_pts_gauss_opt(2,:),bd_pts_gauss_opt(3,:),'.k','markersize',5);
            %[t]=MyCrustOpen(bd_pts');
            %trisurf(t,bd_pts(1,:),bd_pts(2,:),bd_pts(3,:),'facecolor',.5*[1 1 1],'facealpha',.1,'edgecolor','k','edgealpha',.2)
        end
        if strcmp(p.Results.type,'obs')
            plot3(bd_pts_obs_opt(1,:),bd_pts_obs_opt(2,:),bd_pts_obs_opt(3,:),'.','color',.5*[1 1 1],'markersize',5)
        end
        
        %axis image
        grid on
        xlim([min(mu_a(1)-v_a(1,1),mu_b(1)-v_b(1,1)),max(mu_a(1)+v_a(1,1),mu_b(1)+v_b(1,1))])
        ylim([min(mu_a(2)-v_a(2,2),mu_b(2)-v_b(2,2)),max(mu_a(2)+v_a(2,2),mu_b(2)+v_b(2,2))])
        zlim([min(mu_a(3)-v_a(3,3),mu_b(3)-v_b(3,3)),max(mu_a(3)+v_a(3,3),mu_b(3)+v_b(3,3))])

        xlabel('X')
        ylabel('Y')
        zlabel('Z')       

    end
    hold off
end