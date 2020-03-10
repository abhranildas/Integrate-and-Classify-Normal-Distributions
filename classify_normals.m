function results=classify_normals(dist_a,dist_b,varargin)
% Classification between two distributions.
%
% How to use this command:
% See github readme at https://github.com/abhranildas/classify
%
% Credits:
%   Abhranil Das <abhranil.das@utexas.edu>
%	R Calen Walshe
%	Wilson S Geisler
%	Center for Perceptual Systems, University of Texas at Austin
% If you use this code, please cite:
%   A new method to compute classification error
%   https://jov.arvojournals.org/article.aspx?articleid=2750251

% parse inputs
parser = inputParser;
addRequired(parser,'dist_a',@(x) isnumeric(x));
addRequired(parser,'dist_b',@(x) isnumeric(x));
addParameter(parser,'prior_a',0.5, @(x) isnumeric(x) && isscalar(x) && (x > 0) && (x < 1));
addParameter(parser,'vals',eye(2), @(x) isnumeric(x) && ismatrix(x));
addParameter(parser,'custom_bd_coeffs',[]);
addParameter(parser,'custom_bd_fns',[]);
addParameter(parser,'type','params', @(s) strcmp(s,'params') || strcmp(s,'obs'));
addParameter(parser,'bPlot',true, @(x) islogical(x));

parse(parser,dist_a,dist_b,varargin{:});

if strcmp(parser.Results.type,'params')
    mu_a=dist_a(:,1);
    v_a=dist_a(:,2:end);
    mu_b=dist_b(:,1);
    v_b=dist_b(:,2:end);
elseif strcmp(parser.Results.type,'obs')
    mu_a=mean(dist_a)';
    v_a=cov(dist_a);
    mu_b=mean(dist_b)';
    v_b=cov(dist_b);    
end

% if input obs and prior_a is not specified,
if strcmp(parser.Results.type,'obs') && any(strcmp(parser.UsingDefaults,'prior_a'))
    % set prior_a according to input sample sizes
    prior_a=size(dist_a,1)/(size(dist_a,1)+size(dist_b,1));
else
    prior_a=parser.Results.prior_a;
end
prior_b=1-prior_a;

vals=parser.Results.vals;

% check if optimal case
if prior_a==0.5 && isequal(vals,eye(2)) && isempty(parser.Results.custom_bd_coeffs)
    optimal_case=true;
else
    optimal_case=false;
end

dim=length(mu_a); % dimension

if parser.Results.bPlot && dim<=3
    figure; hold on;
end

if ~isempty(parser.Results.custom_bd_fns) % arbitrary boundary functions
    bd_fn_a=parser.Results.custom_bd_fns{1};
    bd_fn_b=parser.Results.custom_bd_fns{2};
else
    if ~isempty(parser.Results.custom_bd_coeffs) % custom boundary coefficients
        bd_coeffs_norm_a=parser.Results.custom_bd_coeffs;
    else
        % compute boundary coefficients
        a2_norm=inv(v_b)-inv(v_a);
        a1_norm=2*(v_a\mu_a-v_b\mu_b);
        a0_norm=mu_b'/v_b*mu_b-mu_a'/v_a*mu_a+log((((vals(1,1)-vals(1,2))*prior_a)/((vals(2,2)-vals(2,1))*prior_b))^2*det(v_b)/det(v_a));
        bd_coeffs_norm_a=struct;
        bd_coeffs_norm_a.a2=a2_norm;
        bd_coeffs_norm_a.a1=a1_norm;
        bd_coeffs_norm_a.a0=a0_norm;
        results.bd_coeffs_norm=bd_coeffs_norm_a;
    end
    % flip boundary sign for b
    bd_coeffs_norm_b=struct;
    bd_coeffs_norm_b.a2=-bd_coeffs_norm_a.a2;
    bd_coeffs_norm_b.a1=-bd_coeffs_norm_a.a1;
    bd_coeffs_norm_b.a0=-bd_coeffs_norm_a.a0;
end

bd_pts_norm=[];
if isequal(mu_a,mu_b) && isequal(v_a,v_b)
    % if the dists are identical:
    acc_norm=0.5;
    acc_norm_a=0.5;
    acc_norm_b=0.5;
else
    % compute accuracy and boundary for each normal, and combine
    if isempty(parser.Results.custom_bd_fns) % if calculated or custom coefficients
        [acc_norm_a,err_norm_a,bd_pts_norm_a]=integrate_normal(mu_a,v_a,'bd_coeffs',bd_coeffs_norm_a,'p_prior',prior_a,'bPlot',parser.Results.bPlot,'plot_color','blue');
        [acc_norm_b,err_norm_b,bd_pts_norm_b]=integrate_normal(mu_b,v_b,'bd_coeffs',bd_coeffs_norm_b,'p_prior',prior_b,'bPlot',parser.Results.bPlot,'plot_color','red');
    else % if arbitrary boundary functions
        [acc_norm_a,err_norm_a,bd_pts_norm_a]=integrate_normal(mu_a,v_a,'bd_fn',bd_fn_a,'p_prior',prior_a,'bPlot',parser.Results.bPlot,'plot_color','blue');
        [acc_norm_b,err_norm_b,bd_pts_norm_b]=integrate_normal(mu_b,v_b,'bd_fn',bd_fn_b,'p_prior',prior_b,'bPlot',parser.Results.bPlot,'plot_color','red');
    end
    acc_norm=prior_a*acc_norm_a+prior_b*acc_norm_b;
    err_norm=prior_a*err_norm_a+prior_b*err_norm_b;
    bd_pts_norm=[bd_pts_norm_a,bd_pts_norm_b];
end

results.bd_pts_norm=bd_pts_norm;
results.errmat_norm=[acc_norm_a, err_norm_a; err_norm_b, acc_norm_b];
results.err_norm=err_norm;

% d'
if optimal_case
    d_norm=-2*norminv(err_norm);
    results.d_norm=d_norm;
end

if ~isequal(vals,eye(2)) % if outcome values are supplied
    results.outcome_vals_norm=[vals(1,1)*prior_a*acc_norm_a, vals(1,2)*prior_a*(1-acc_norm_a);...
        vals(2,1)*prior_b*(1-acc_norm_b), vals(2,2)*prior_b*acc_norm_b];
    results.ex_val_norm=(vals(1,1)-vals(1,2))*prior_a*acc_norm_a + ...
        (vals(2,2)-vals(2,1))*prior_b*acc_norm_b + vals(1,2)*prior_a + vals(2,1)*prior_b;
end
    
% if input is observations,
if strcmp(parser.Results.type,'obs')
    % if custom boundary is provided,
    if ~isempty(parser.Results.custom_bd_coeffs)
        % compute accuracy and expected value with custom boundary
        [acc_obs,acc_obs_a,acc_obs_b,outcome_counts_obs]=val_obs(dist_a,dist_b,bd_coeffs_norm_a,eye(size(dist_a,2)));
        [ex_val_obs,~,~,outcome_vals_obs]=val_obs(dist_a,dist_b,bd_coeffs_norm_a,vals);
    else
        %acc_obs_start=accuracy_obs(dist_a,dist_b,opt_bd_coeffs_norm);
        % find boundary that optimizes expected value / accuracy
        fun = @(x)-val_obs_flat(dim,x,dist_a,dist_b,vals);
        x=fminsearch(fun,[a2_norm(:); a1_norm(:); a0_norm],optimset('Display','iter'));

        a2_obs=reshape(x(1:dim^2),[dim dim])';
        a1_obs=x(dim^2+1:dim^2+dim);
        a0_obs=x(end);

        bd_coeffs_obs_a=struct;
        bd_coeffs_obs_a.a2=a2_obs;
        bd_coeffs_obs_a.a1=a1_obs;
        bd_coeffs_obs_a.a0=a0_obs;
        results.bd_coeffs_obs=bd_coeffs_obs_a;
        
        % flip boundary sign for b
        bd_coeffs_obs_b.a2=-bd_coeffs_obs_a.a2;
        bd_coeffs_obs_b.a1=-bd_coeffs_obs_a.a1;
        bd_coeffs_obs_b.a0=-bd_coeffs_obs_a.a0;
        
        if dim<=3
            % boundary points
            [~,~,bd_pts_obs_a]=integrate_normal(mu_a,v_a,'bd_coeffs',bd_coeffs_obs_a,'bPlot',false);
            [~,~,bd_pts_obs_b]=integrate_normal(mu_b,v_b,'bd_coeffs',bd_coeffs_obs_b,'bPlot',false);
            bd_pts_obs=[bd_pts_obs_a,bd_pts_obs_b];
            results.bd_pts_obs=bd_pts_obs;
        end
        
        % accuracy with optimal data boundary
        [acc_obs,acc_obs_a,acc_obs_b,outcome_counts_obs]=val_obs(dist_a,dist_b,bd_coeffs_obs_a,eye(2));
        
        % expected value with optimal data boundary
        [ex_val_obs,~,~,outcome_vals_obs]=val_obs(dist_a,dist_b,bd_coeffs_obs_a,vals);
        
        if optimal_case && size(dist_a,1)==size(dist_b,1) % if both category frequencies same,
            d_obs=2*norminv(acc_obs); % obs accuracy can be used to compute obs d'
            results.d_obs=d_obs;
        end
    end
    results.errmat_obs=[acc_obs_a, 1-acc_obs_a; 1-acc_obs_b, acc_obs_b];
    results.err_obs=1-acc_obs;
    results.outcome_counts_obs=outcome_counts_obs;
    if ~isequal(vals,eye(2)) % if outcome values are supplied
        results.ex_val_obs=ex_val_obs;
        results.outcome_vals_obs=outcome_vals_obs;
    end
end

% approximate d' and Chernoff bound
if optimal_case
    % approximate d' (considering each vcov = avg)
    vcov_av=(v_a+v_b)/2;
    results.d_norm_aprx=sqrt((mu_a-mu_b)'/(vcov_av)*(mu_a-mu_b));    
    % Chernoff bound for d'
    [~,k]=fminbnd(@(b)chernoff_bound(b,mu_a,v_a,mu_b,v_b,0.5),0,1);
    err_opt_min_chernoff=10^k;
    d_norm_min=-2*norminv(err_opt_min_chernoff);
    results.d_norm_min=d_norm_min;
end
% Chernoff bound for error
[~,k]=fminbnd(@(b)chernoff_bound(b,mu_a,v_a,mu_b,v_b,prior_a),0,1);
results.log_err_norm_max=k;

%% Plot data
if parser.Results.bPlot && dim<=3
    hold on;
    % plot title
    if strcmp(parser.Results.type,'params')
        title(sprintf("error = %g",1-acc_norm(1)))
    else
        title(sprintf("error = %g",1-acc_obs(1)))
    end
    
    if dim==1
        % plot data points
        if strcmp(parser.Results.type,'obs')
            [heights_a,edges_a]=histcounts(dist_a,'BinMethod','scott','normalization','pdf');
            histogram('BinCounts',heights_a*prior_a,'BinEdges',edges_a,'facecolor','blue','facealpha',0.3,'edgecolor','blue','edgealpha',0.3);
            [heights_b,edges_b]=histcounts(dist_b,'BinMethod','scott','normalization','pdf');
            histogram('BinCounts',heights_b*prior_b,'BinEdges',edges_b,'facecolor','red','facealpha',0.3,'edgecolor','red','edgealpha',0.3);
        end
        % plot data boundary
        if strcmp(parser.Results.type,'obs')&& isempty(parser.Results.custom_bd_coeffs)
            for x=bd_pts_obs
                line([x x],ylim,'color',.5*[1 1 1],'linewidth',1)
            end
        end        
    elseif dim==2
        % plot data points
        if strcmp(parser.Results.type,'obs')
            scatter(dist_a(:,1),dist_a(:,2),4,'o','markerfacecolor','blue');
            scatter(dist_b(:,1),dist_b(:,2),4,'o','markerfacecolor','red');
        end
        % plot data boundary
        if strcmp(parser.Results.type,'obs') && isempty(parser.Results.custom_bd_coeffs)
            plot(bd_pts_obs(1,:),bd_pts_obs(2,:),'.','color',.5*[1 1 1])
        end
    elseif dim==3
        % plot data points
        if strcmp(parser.Results.type,'obs')
            plot3(dist_a(:,1),dist_a(:,2),dist_a(:,3),'.','color','blue','markersize',5);
            plot3(dist_b(:,1),dist_b(:,2),dist_b(:,3),'.','color','red','markersize',5);
        end
        % plot data boundary
        if strcmp(parser.Results.type,'obs') && isempty(parser.Results.custom_bd_coeffs)
            plot3(bd_pts_obs(1,:),bd_pts_obs(2,:),bd_pts_obs(3,:),'.','color',.5*[1 1 1],'markersize',5)
        end
        xlim([min(mu_a(1)-v_a(1,1),mu_b(1)-v_b(1,1)),max(mu_a(1)+v_a(1,1),mu_b(1)+v_b(1,1))])
        ylim([min(mu_a(2)-v_a(2,2),mu_b(2)-v_b(2,2)),max(mu_a(2)+v_a(2,2),mu_b(2)+v_b(2,2))])
        zlim([min(mu_a(3)-v_a(3,3),mu_b(3)-v_b(3,3)),max(mu_a(3)+v_a(3,3),mu_b(3)+v_b(3,3))])
    end
    hold off
end