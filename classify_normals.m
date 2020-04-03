function results=classify_normals(dist_1,dist_2,varargin)
% Classification between two distributions.
%
% How to use this command:
% See github readme at https://github.com/abhranildas/classify
%
% Author:
%   Abhranil Das <abhranil.das@utexas.edu>
%	Center for Perceptual Systems, University of Texas at Austin
% If you use this code, please cite:
%   A new method to compute classification error
%   https://jov.arvojournals.org/article.aspx?articleid=2750251

% parse inputs
parser = inputParser;
addRequired(parser,'dist_1',@isnumeric);
addRequired(parser,'dist_2',@isnumeric);
addParameter(parser,'prior_1',0.5, @(x) isnumeric(x) && isscalar(x) && (x > 0) && (x < 1));
addParameter(parser,'vals',eye(2), @(x) isnumeric(x) && ismatrix(x));
addParameter(parser,'reg',[]);
addParameter(parser,'reg_type','quad');
addParameter(parser,'type','norm', @(s) strcmp(s,'norm') || strcmp(s,'samp'));
addParameter(parser,'n_rays',1e4,@isnumeric);
addParameter(parser,'estimate',[],@(x) (isnumeric(x)&&isscalar(x))||strcmpi(x,'tail'));
addParameter(parser,'bPlot',true, @islogical);

parse(parser,dist_1,dist_2,varargin{:});
reg=parser.Results.reg;
reg_type=parser.Results.reg_type;
estimate=parser.Results.estimate;

if strcmp(parser.Results.type,'norm')
    mu_1=dist_1(:,1);
    v_1=dist_1(:,2:end);
    mu_2=dist_2(:,1);
    v_2=dist_2(:,2:end);
elseif strcmp(parser.Results.type,'samp')
    mu_1=mean(dist_1)';
    v_1=cov(dist_1);
    mu_2=mean(dist_2)';
    v_2=cov(dist_2);
end

% if input samp and prior is not specified,
if strcmp(parser.Results.type,'samp') && any(strcmp(parser.UsingDefaults,'prior_1'))
    % set prior according to input sample sizes
    priors(1)=size(dist_1,1)/(size(dist_1,1)+size(dist_2,1));
else
    priors(1)=parser.Results.prior_1;
end
priors(2)=1-priors(1);

vals=parser.Results.vals;
n_rays=parser.Results.n_rays;
bPlot=parser.Results.bPlot;

% check if optimal case
if priors(1)==0.5 && isequal(vals,eye(2)) && isempty(reg)
    optimal_case=true;
else
    optimal_case=false;
end

dim=length(mu_1); % dimension

colors=colororder;
if bPlot && dim<=3, figure, hold on, end

norm_bd_pts=[];
if isequal(mu_1,mu_2) && isequal(v_1,v_2)
    % if the dists are identical:
    norm_acc_1=0.5;
    norm_acc_2=0.5;
else
    % compute accuracy and boundary for each normal, and combine
    if strcmp(reg_type,'quad') % optimal or custom quad coefficients
        if ~isempty(reg) % custom quad coefficients
            norm_reg_quad_1=reg;
        else
            % compute optimal boundary coefficients
            norm_reg_quad_1=opt_reg_quad([mu_1,v_1],[mu_2,v_2],'vals',vals,'prior_1',priors(1));
            results.norm_reg_quad=norm_reg_quad_1;
        end
        % flip boundary sign for 2nd normal
        norm_reg_quad_2=struct;
        norm_reg_quad_2.a2=-norm_reg_quad_1.a2;
        norm_reg_quad_2.a1=-norm_reg_quad_1.a1;
        norm_reg_quad_2.a0=-norm_reg_quad_1.a0;
        
        [norm_acc_1,norm_err_1,norm_bd_pts_a]=integrate_normal(mu_1,v_1,norm_reg_quad_1,'reg_type','quad','prior',priors(1),'n_rays',n_rays,'estimate',estimate,'bPlot',bPlot,'plot_color',colors(1,:));
        if bPlot && dim<=3, hold on, end
        [norm_acc_2,norm_err_b,norm_bd_pts_b]=integrate_normal(mu_2,v_2,norm_reg_quad_2,'reg_type','quad','prior',priors(2),'n_rays',n_rays,'estimate',estimate,'bPlot',bPlot,'plot_color',colors(2,:));
    elseif strcmp(reg_type,'ray_scan') % ray-scanned region functions
        [norm_acc_1,norm_err_1,norm_bd_pts_a]=integrate_normal(mu_1,v_1,reg{1},'reg_type','ray_scan','prior',priors(1),'n_rays',n_rays,'bPlot',bPlot,'plot_color',colors(1,:));
        if bPlot && dim<=3, hold on, end
        [norm_acc_2,norm_err_b,norm_bd_pts_b]=integrate_normal(mu_2,v_2,reg{2},'reg_type','ray_scan','prior',priors(2),'n_rays',n_rays,'bPlot',bPlot,'plot_color',colors(2,:));
    end
    norm_err=priors(1)*norm_err_1+priors(2)*norm_err_b;
    norm_bd_pts=uniquetol([norm_bd_pts_a,norm_bd_pts_b]',1e-12,'Byrows',true,'Datascale',1)'; % trim to unique boundary points
end

if ~isempty(norm_bd_pts)
    results.norm_bd_pts=norm_bd_pts;
end
results.norm_err_mat=[norm_acc_1, norm_err_1; norm_err_b, norm_acc_2];
results.norm_err=norm_err;

% d'
if optimal_case
    norm_d=-2*norminv(norm_err);
    results.norm_d=norm_d;
end

% Mahalanobis d' and Chernoff bound
if optimal_case
    % approximate d' (considering each vcov = avg)
    vcov_av=(v_1+v_2)/2;
    results.norm_maha_d=sqrt((mu_1-mu_2)'/(vcov_av)*(mu_1-mu_2));
    % Chernoff bound for d'
    [~,k]=fminbnd(@(b)chernoff_bound(b,mu_1,v_1,mu_2,v_2,0.5),0,1);
    err_opt_min_chernoff=10^k;
    norm_min_d=-2*norminv(err_opt_min_chernoff);
    results.norm_min_d=norm_min_d;
end
% Chernoff bound for error
[~,k]=fminbnd(@(b)chernoff_bound(b,mu_1,v_1,mu_2,v_2,priors(1)),0,1);
results.norm_log_max_err=k;

if ~isequal(vals,eye(2)) % if outcome values are supplied
    results.norm_val_mat=results.norm_err_mat.*vals; % conditional expected values
    results.norm_val=sum(sum(results.norm_val_mat.*priors'));
end

% if sample input,
if strcmp(parser.Results.type,'samp')
    if strcmp(reg_type,'quad') % if default or custom quad boundary,
        % compute outcome counts and error
        [~,samp_count_mat]=value_samp(dist_1,dist_2,norm_reg_quad_1,ones(2));
        samp_err=value_samp(dist_1,dist_2,norm_reg_quad_1,~eye(2));
        results.samp_count_mat=samp_count_mat;
        results.samp_err_mat=samp_count_mat./sum(samp_count_mat,2);
        results.samp_err=samp_err;
        
        if optimal_case
            results.samp_d=-2*norminv(samp_err); % samp error can be used to compute samp d'
        elseif ~isequal(vals,eye(2)) % if outcome values are supplied
            [samp_ex_val,samp_val_mat]=value_samp(dist_1,dist_2,norm_reg_quad_1,vals);
            results.samp_val_mat=samp_val_mat;
            results.samp_ex_val=samp_ex_val;
        end            

        if isempty(reg) % if no custom boundary
            % find quad boundary that optimizes expected value / accuracy
            x=fminsearch(@(x)-value_samp_flat(dim,x,dist_1,dist_2,vals),[norm_reg_quad_1.a2(:); norm_reg_quad_1.a1(:); norm_reg_quad_1.a0],optimset('Display','iter'));
            
            a2_samp=reshape(x(1:dim^2),[dim dim])';
            a1_samp=x(dim^2+1:dim^2+dim);
            a0_samp=x(end);
            
            samp_opt_reg_quad_1=struct;
            samp_opt_reg_quad_1.a2=a2_samp;
            samp_opt_reg_quad_1.a1=a1_samp;
            samp_opt_reg_quad_1.a0=a0_samp;
            results.samp_opt_reg_quad=samp_opt_reg_quad_1;
            
            % flip boundary sign for b
            samp_opt_reg_quad_2.a2=-samp_opt_reg_quad_1.a2;
            samp_opt_reg_quad_2.a1=-samp_opt_reg_quad_1.a1;
            samp_opt_reg_quad_2.a0=-samp_opt_reg_quad_1.a0;
            
            if dim<=3
                % boundary points
                [~,~,samp_opt_bd_pts_1]=integrate_normal(mu_1,v_1,samp_opt_reg_quad_1,'n_rays',n_rays,'bPlot',false);
                [~,~,samp_opt_bd_pts_2]=integrate_normal(mu_2,v_2,samp_opt_reg_quad_2,'n_rays',n_rays,'bPlot',false);
                samp_opt_bd_pts=[samp_opt_bd_pts_1,samp_opt_bd_pts_2];
                results.samp_opt_bd_pts=samp_opt_bd_pts;
            end
            
            % compute outcome counts and error with optimized data boundary
            [~,samp_opt_count_mat]=value_samp(dist_1,dist_2,samp_opt_reg_quad_1,ones(2));
            samp_opt_err=value_samp(dist_1,dist_2,samp_opt_reg_quad_1,~eye(2));
            results.samp_opt_count_mat=samp_opt_count_mat;
            results.samp_opt_err_mat=samp_opt_count_mat./sum(samp_opt_count_mat,2);
            results.samp_opt_err=samp_opt_err;
            
            if optimal_case
                results.samp_opt_d=-2*norminv(samp_opt_err); % samp error can be used to compute samp d'
            elseif ~isequal(vals,eye(2)) % if outcome values are supplied
                [samp_opt_ex_val,samp_opt_val_mat]=value_samp(dist_1,dist_2,samp_opt_reg_quad_1,vals);
                results.samp_opt_val_mat=samp_opt_val_mat;
                results.samp_opt_ex_val=samp_opt_ex_val;
            end
        end
    end
end

%% Plot samples
if bPlot && dim<=3
    hold on    
    if strcmp(parser.Results.type,'norm')
        title(sprintf("error = %g",norm_err)) % plot title
    elseif strcmp(parser.Results.type,'samp')
        if ~isempty(reg) % if custom boundary
            title(sprintf("error = %g / %g",[norm_err,samp_err])) % plot title
            % don't plot sample boundary
            plot_sample(dist_1,[],priors(1),colors(1,:))
            plot_sample(dist_2,[],priors(2),colors(2,:))
        else
            title(sprintf("error = %g / %g / %g",[norm_err,samp_err,samp_opt_err])) % plot title
            plot_sample(dist_1,samp_opt_bd_pts_1,priors(1),colors(1,:))
            plot_sample(dist_2,samp_opt_bd_pts_2,priors(2),colors(2,:))
            
            % boundary legends
            if dim<=2
                norm_bd_marker=line(nan, nan,'color','k','linewidth',1);
                samp_bd_marker=line(nan, nan, 'color', .5*[1 1 1],'linewidth',1);
                
            else
                norm_bd_marker=line(nan, nan,'color','k','linewidth',1.5,'linestyle',':');
                samp_bd_marker=line(nan, nan, 'color', .5*[1 1 1],'linewidth',1.5,'linestyle',':');
            end
            legend([norm_bd_marker,samp_bd_marker],'normal boundary','sample boundary','edgecolor','none')
        end
    end
    hold off
end