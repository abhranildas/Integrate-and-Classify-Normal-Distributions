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

%% parse inputs
parser = inputParser;
addRequired(parser,'dist_1',@isnumeric);
addRequired(parser,'dist_2',@isnumeric);
addParameter(parser,'prior_1',0.5, @(x) isnumeric(x) && isscalar(x) && (x > 0) && (x < 1));
addParameter(parser,'vals',eye(2), @(x) isnumeric(x) && ismatrix(x));
addParameter(parser,'reg',[]);
addParameter(parser,'reg_type','quad');
addParameter(parser,'cheb_reg_span',3);
addParameter(parser,'func_crossings',100);
addParameter(parser,'type','norm', @(s) strcmp(s,'norm') || strcmp(s,'samp'));
addParameter(parser,'opt_samp',true, @islogical);
addParameter(parser,'AbsTol',1e-10);
addParameter(parser,'RelTol',1e-2);
addParameter(parser,'n_samp_bd_pts',1e4);
addParameter(parser,'bPlot',true, @islogical);

parse(parser,dist_1,dist_2,varargin{:});
reg=parser.Results.reg;
reg_type=parser.Results.reg_type;
cheb_reg_span=parser.Results.cheb_reg_span;
func_crossings=parser.Results.func_crossings;
AbsTol=parser.Results.AbsTol;
RelTol=parser.Results.RelTol;
vals=parser.Results.vals;
n_samp_bd_pts=parser.Results.n_samp_bd_pts;
bPlot=parser.Results.bPlot;

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

% check if optimal case
if priors(1)==0.5 && isequal(vals,eye(2)) && isempty(reg)
    optimal_case=true;
else
    optimal_case=false;
end

dim=length(mu_1); % dimension

colors=colororder;
if bPlot, figure, hold on, end

%% normal inputs

norm_bd_pts=[];
if isequal(mu_1,mu_2) && isequal(v_1,v_2)
    % if the dists are identical:
    norm_err_1=0.5; norm_acc_1=0.5;
    norm_err_2=0.5; norm_acc_2=0.5;
else
    % compute accuracy and boundary for each normal, and combine
    if isempty(reg) % optimal region
        % compute optimal boundary coefficients
        reg=opt_reg_quad([mu_1,v_1],[mu_2,v_2],'vals',vals,'prior_1',priors(1));
        results.norm_reg_quad=reg;
    end
    if strcmp(reg_type,'quad') % custom quad coefficients
        % flip boundary sign for 2nd normal
        norm_reg_2=structfun(@uminus,reg,'un',0);
        
        [norm_acc_1,norm_err_1,norm_bd_pts_1]=integrate_normal(mu_1,v_1,reg,'prior',priors(1),'AbsTol',AbsTol,'RelTol',RelTol,'bPlot',bPlot,'plot_color',colors(1,:));
        if bPlot && dim<=3, hold on, end
        [norm_acc_2,norm_err_2,norm_bd_pts_2]=integrate_normal(mu_2,v_2,norm_reg_2,'prior',priors(2),'AbsTol',AbsTol,'RelTol',RelTol,'bPlot',bPlot,'plot_color',colors(2,:));
        
    elseif strcmp(reg_type,'ray_scan') % ray-scanned region functions
        [norm_acc_1,norm_err_1,norm_bd_pts_1]=integrate_normal(mu_1,v_1,reg{1},'reg_type','ray_scan','prior',priors(1),'AbsTol',AbsTol,'RelTol',RelTol,'bPlot',bPlot,'plot_color',colors(1,:));
        if bPlot && dim<=3, hold on, end
        [norm_acc_2,norm_err_2,norm_bd_pts_2]=integrate_normal(mu_2,v_2,reg{2},'reg_type','ray_scan','prior',priors(2),'AbsTol',AbsTol,'RelTol',RelTol,'bPlot',bPlot,'plot_color',colors(2,:));
        
    elseif strcmp(reg_type,'cheb') % cheb region function
        [norm_acc_1,norm_err_1,norm_bd_pts_1]=integrate_normal(mu_1,v_1,reg,'reg_type','cheb','cheb_reg_span',cheb_reg_span,'func_crossings',func_crossings,'prior',priors(1),'AbsTol',AbsTol,'RelTol',RelTol,'bPlot',bPlot,'plot_color',colors(1,:));
        if bPlot && dim<=3, hold on, end
        [norm_err_2,norm_acc_2,norm_bd_pts_2]=integrate_normal(mu_2,v_2,reg,'reg_type','cheb','cheb_reg_span',cheb_reg_span,'func_crossings',func_crossings,'prior',priors(2),'AbsTol',AbsTol,'RelTol',RelTol,'bPlot',bPlot,'plot_color',colors(2,:));
    end
    norm_bd_pts=uniquetol([norm_bd_pts_1,norm_bd_pts_2]',1e-12,'Byrows',true,'Datascale',1)'; % trim to unique boundary points
end
norm_err=priors(1)*norm_err_1+priors(2)*norm_err_2;

if ~isempty(norm_bd_pts)
    results.norm_bd_pts=norm_bd_pts;
end

results.norm_err_mat=[norm_acc_1, norm_err_1; norm_err_2, norm_acc_2];
results.norm_err=norm_err;

% d'
if optimal_case
    norm_dprime=-2*norminv(norm_err);
    results.norm_dprime=norm_dprime;
end

% Mahalanobis d'
if optimal_case
    % approximate d' (considering each vcov = avg)
    vcov_av=(v_1+v_2)/2;
    results.norm_maha_dprime=sqrt((mu_1-mu_2)'/(vcov_av)*(mu_1-mu_2));
end

if ~isequal(vals,eye(2)) % if outcome values are supplied
    results.norm_val_mat=results.norm_err_mat.*vals; % conditional expected values
    results.norm_val=sum(sum(results.norm_val_mat.*priors'));
end

%% sample inputs
if strcmp(parser.Results.type,'samp')
    
    % compute outcome counts and error
    [~,samp_count_mat]=samp_value(dist_1,dist_2,reg,'reg_type',reg_type,'vals',ones(2));
    samp_err=samp_value(dist_1,dist_2,reg,'reg_type',reg_type,'vals',1-eye(2));
    
    if optimal_case
        results.samp_dprime=-2*norminv(samp_err); % samp error can be used to compute samp d'
    elseif ~isequal(vals,eye(2)) % if outcome values are supplied
        [samp_ex_val,samp_val_mat]=samp_value(dist_1,dist_2,reg,'reg_type',reg_type,'vals',vals);
        results.samp_val_mat=samp_val_mat;
        results.samp_ex_val=samp_ex_val;
    end
    
    results.samp_count_mat=samp_count_mat;
    results.samp_err_mat=samp_count_mat./sum(samp_count_mat,2);
    results.samp_err=samp_err;
    
    % log posterior ratios for >3D
    lpr_1=normal_lpr(dist_1,mu_1',v_1,mu_2',v_2,priors(1));
    lpr_2=normal_lpr(dist_2,mu_1',v_1,mu_2',v_2,priors(1));
    results.samp_lpr={lpr_1,lpr_2};
    
    %% sample-optimal boundary
    if isempty(parser.Results.reg) && parser.Results.opt_samp % if default boundary,
        try
            % find quad boundary that optimizes expected value / accuracy
            x=fminsearch(@(x) -samp_value_flat(dim,x,dist_1,dist_2,vals),[reg.a2(:); reg.a1(:); reg.a0],optimset('Display','iter','TolFun',1e-6));
            
            samp_reg_1.a2=reshape(x(1:dim^2),[dim dim])';
            samp_reg_1.a1=x(dim^2+1:dim^2+dim);
            samp_reg_1.a0=x(end);
            results.samp_opt_reg_quad=samp_reg_1;
            
            % flip boundary sign for b
            samp_reg_2=structfun(@uminus,samp_reg_1,'un',0);
            
            if dim<=3
                % boundary points
                [~,samp_bd_pts_1]=prob_bd_angle(mu_1,v_1,samp_reg_1,'n_bd_pts',n_samp_bd_pts);
                [~,samp_bd_pts_2]=prob_bd_angle(mu_2,v_2,samp_reg_2,'n_bd_pts',n_samp_bd_pts);
                samp_bd_pts=uniquetol([samp_bd_pts_1,samp_bd_pts_2]',1e-12,'Byrows',true,'Datascale',1)'; % trim to unique boundary points
                results.samp_opt_bd_pts=samp_bd_pts;
            end
            
            samp_opt_flag=true;
            
        catch mException
            %warning(getReport(mException,'extended','hyperlinks','on'))
            warning('Cannot optimize full quadratic classification boundary for the sample (usually due to system limits). Finding best log posterior ratio classification criterion instead.')
            samp_opt_flag=false;
            
            if dim>3
                % find optimal lpr
                samp_opt_lpr=-fminbnd(@(x) -samp_value_flat(1,[0;1;x],lpr_1,lpr_2,vals),min(lpr_2),max(lpr_1),optimset('Display','iter'));
                %x=fminsearch(@(x) -samp_value_flat(1,[0;1;x],lpr_1,lpr_2,vals),0,optimset('TolFun',1e-6,'Display','iter'));
                results.samp_opt_lpr=samp_opt_lpr;
                
                samp_reg_1.a2=0;
                samp_reg_1.a1=1;
                samp_reg_1.a0=-samp_opt_lpr;
                dist_1=lpr_1;
                dist_2=lpr_2;
            end
        end
        
        % compute outcome counts and error with optimized sample boundary
        [~,samp_opt_count_mat]=samp_value(dist_1,dist_2,samp_reg_1,'vals',ones(2));
        samp_opt_err=samp_value(dist_1,dist_2,samp_reg_1,'vals',1-eye(2));
        results.samp_opt_count_mat=samp_opt_count_mat;
        results.samp_opt_err_mat=samp_opt_count_mat./sum(samp_opt_count_mat,2);
        results.samp_opt_err=samp_opt_err;
        
        if optimal_case
            results.samp_opt_dprime=-2*norminv(samp_opt_err); % samp error can be used to compute samp d'
        elseif ~isequal(vals,eye(2)) % if outcome values are supplied
            [samp_opt_ex_val,samp_opt_val_mat]=samp_value(dist_1,dist_2,samp_reg_1,'vals',vals);
            results.samp_opt_val_mat=samp_opt_val_mat;
            results.samp_opt_ex_val=samp_opt_ex_val;
        end
    end
    
    
end

%% Plot samples
if bPlot
    hold on
    
    % plot log posterior ratio distributions
    if dim>3
        lpr_quad=opt_reg_quad([mu_1,v_1],[mu_2,v_2],'prior_1',priors(1));
        
        if nnz(lpr_quad.a2) % if the quadratic term exists, i.e. unequal vcovs
            % lpr distributions are generalized chi-squared
            [lambda_1,m_1,delta_1,c_1]=norm_quad_to_gx2_params(mu_1,v_1,lpr_quad);
            [mu_llr_1,v_llr_1]=gx2stat(lambda_1,m_1,delta_1,c_1); % mean and variance of llr
            x_1=linspace(mu_llr_1-5*sqrt(v_llr_1),mu_llr_1+5*sqrt(v_llr_1),1e3);
            y_1=priors(1)*arrayfun(@(x) gx2pdf(x,lambda_1,m_1,delta_1,c_1),x_1);
            area(x_1,y_1,'facecolor',colors(1,:),'facealpha',0.4,'edgecolor',colors(1,:),'edgealpha',0.5,'linewidth',1)
            
            [lambda_2,m_2,delta_2,c_2]=norm_quad_to_gx2_params(mu_2,v_2,lpr_quad);
            [mu_llr_2,v_llr_2]=gx2stat(lambda_2,m_2,delta_2,c_2); % mean and variance of llr
            x_2=linspace(mu_llr_2-5*sqrt(v_llr_2),mu_llr_2+5*sqrt(v_llr_2),1e3);
            y_2=priors(2)*arrayfun(@(x) gx2pdf(x,lambda_2,m_2,delta_2,c_2),x_2);
            area(x_2,y_2,'facecolor',colors(2,:),'facealpha',0.4,'edgecolor',colors(2,:),'edgealpha',0.5,'linewidth',1)
            
        else % lpr distributions are normal
            mu_llr_1=lpr_quad.a1'*mu_1+lpr_quad.a0;
            v_llr_1=lpr_quad.a1'*v_1*lpr_quad.a1;
            plot_normal(mu_llr_1,v_llr_1,priors(1),colors(1,:)); hold on
            
            mu_llr_2=lpr_quad.a1'*mu_2+lpr_quad.a0;
            v_llr_2=lpr_quad.a1'*v_2*lpr_quad.a1;
            plot_normal(mu_llr_2,v_llr_2,priors(2),colors(2,:)); hold on
        end
        % plot lpr boundary
        xline(-log((vals(1,1)-vals(1,2))/(vals(2,2)-vals(2,1))),'color','k','linewidth',1);
        xlabel('log posterior ratio: $\log \frac{p(x \epsilon N_1)}{p(x \epsilon N_2)}$','interpreter','latex','fontsize',15)
    end
    
    if strcmp(parser.Results.type,'norm')
        title(sprintf("error = %g",norm_err)) % plot title
    elseif strcmp(parser.Results.type,'samp')
        
        if dim <=3
            plot_sample(dist_1,priors(1),colors(1,:));            
            plot_sample(dist_2,priors(2),colors(2,:));            
        elseif dim>3
            % plot sample log posterior ratios
            plot_sample(lpr_1,priors(1),colors(1,:))
            plot_sample(lpr_2,priors(2),colors(2,:))
        end
        if ~isempty(parser.Results.reg) % if custom boundary
            title(sprintf("error = %g / %g",[norm_err,samp_err])) % plot title
            % don't plot sample boundary
        else
            if parser.Results.opt_samp
                title(sprintf("error = %g / %g / %g",[norm_err,samp_err,samp_opt_err])) % plot title
                if dim <=3
                    plot_boundary(samp_reg_1,dim,'reg_type','quad','orig',mu_1,'v',v_1,'plot_color',.5*[1 1 1]);
                    plot_boundary(samp_reg_1,dim,'reg_type','quad','orig',mu_2,'v',v_2,'plot_color',.5*[1 1 1]);
                elseif dim>3 && ~samp_opt_flag
                    xline(samp_opt_lpr,'color',.5*[1 1 1]);
                end
            else
                title(sprintf("error = %g / %g",[norm_err,samp_err])) % plot title
            end
        end
        % boundary legends
        if dim==3
            norm_bd_marker=line(nan, nan,'color','k','linewidth',1.5,'linestyle',':');
            samp_bd_marker=line(nan, nan, 'color', .5*[1 1 1],'linewidth',1.5,'linestyle',':');
        else
            norm_bd_marker=line(nan, nan,'color','k','linewidth',1);
            samp_bd_marker=line(nan, nan, 'color', .5*[1 1 1],'linewidth',1);
        end
        if dim>3 && samp_opt_flag
            legend(norm_bd_marker,'normal boundary','edgecolor','none')
        else
            legend([norm_bd_marker,samp_bd_marker],'normal boundary','sample boundary','edgecolor','none')
        end
    end
    hold off
end