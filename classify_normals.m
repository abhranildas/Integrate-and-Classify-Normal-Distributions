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
parser=inputParser;
parser.KeepUnmatched=true;
addRequired(parser,'dist_1',@isnumeric);
addRequired(parser,'dist_2',@isnumeric);
addParameter(parser,'prior_1',0.5, @(x) isnumeric(x) && isscalar(x) && (x > 0) && (x < 1));
addParameter(parser,'vals',eye(2), @(x) isnumeric(x) && ismatrix(x));
addParameter(parser,'reg',[]);
addParameter(parser,'reg_type','quad');
addParameter(parser,'method','ray');
addParameter(parser,'fun_span',3);
addParameter(parser,'fun_resol',100);
addParameter(parser,'type','norm', @(s) strcmpi(s,'norm') || strcmpi(s,'samp'));
addParameter(parser,'samp_opt',true, @islogical);
addParameter(parser,'AbsTol',1e-10);
addParameter(parser,'RelTol',1e-2);
addParameter(parser,'n_samp_bd_pts',1e4);
addParameter(parser,'bPlot',true,@islogical);

parse(parser,dist_1,dist_2,varargin{:});
reg=parser.Results.reg;
reg_type=parser.Results.reg_type;
vals=parser.Results.vals;
n_samp_bd_pts=parser.Results.n_samp_bd_pts;
bPlot=parser.Results.bPlot;

if strcmpi(parser.Results.type,'norm')
    mu_1=dist_1(:,1);
    v_1=dist_1(:,2:end);
    mu_2=dist_2(:,1);
    v_2=dist_2(:,2:end);
elseif strcmpi(parser.Results.type,'samp')
    mu_1=mean(dist_1)';
    v_1=cov(dist_1);
    mu_2=mean(dist_2)';
    v_2=cov(dist_2);
end

% if input samp and prior is not specified,
if strcmpi(parser.Results.type,'samp') && any(strcmpi(parser.UsingDefaults,'prior_1'))
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
        reg=opt_class_quad([mu_1,v_1],[mu_2,v_2],'vals',vals,'prior_1',priors(1));
        results.norm_bd=reg;
    end
    
    if strcmpi(reg_type,'ray_scan') % ray-scanned region functions
        [norm_acc_1,norm_err_1,norm_bd_pts_1]=integrate_normal(mu_1,v_1,reg,'prior',priors(1),'bPlot',bPlot*2,'plot_color',colors(1,:),varargin{:});
        if bPlot && dim<=3, hold on, end
        reg_inv=@(n,mu,v) invert_reg(reg,n,'mu',mu,'v',v);
        [norm_acc_2,norm_err_2,norm_bd_pts_2]=integrate_normal(mu_2,v_2,reg_inv,'prior',priors(2),'bPlot',bPlot*2,'plot_color',colors(2,:),varargin{:});
    else
        if strcmpi(reg_type,'quad') % custom quad coefficients
            % flip boundary sign for 2nd normal
            %norm_reg_2=structfun(@uminus,reg,'un',0);
            [norm_acc_1,norm_err_1,norm_bd_pts_1]=integrate_normal(mu_1,v_1,reg,'prior',priors(1),'plot_color',colors(1,:),varargin{:});
            if bPlot, hold on, end
            [norm_err_2,norm_acc_2,norm_bd_pts_2]=integrate_normal(mu_2,v_2,reg,'prior',priors(2),'plot_color',colors(2,:),varargin{:});
        elseif strcmpi(reg_type,'fun') % region defined by a function
            [norm_acc_1,norm_err_1,norm_bd_pts_1]=integrate_normal(mu_1,v_1,reg,'prior',priors(1),'plot_color',colors(1,:),varargin{:});
            if bPlot && dim<=3, hold on, end
            [norm_err_2,norm_acc_2,norm_bd_pts_2]=integrate_normal(mu_2,v_2,reg,'prior',priors(2),'plot_color',colors(2,:),varargin{:});
        end
        % plot boundary
        if bPlot, hold on, plot_boundary(reg,dim,'reg_type',reg_type), end
    end
    norm_bd_pts=uniquetol([norm_bd_pts_1,norm_bd_pts_2]',1e-12,'Byrows',true,'Datascale',1)'; % trim to unique boundary points
end
%norm_err=priors(1)*norm_err_1+priors(2)*norm_err_2;

if ~isempty(norm_bd_pts)
    results.norm_bd_pts=norm_bd_pts;
end

norm_errmat=[[norm_acc_1, norm_err_1]*priors(1); [norm_err_2, norm_acc_2]*priors(2)];
results.norm_errmat=norm_errmat;
%results.norm_errs=[norm_acc_1, norm_err_1; norm_err_2, norm_acc_2];
norm_err=sum(norm_errmat(~eye(2))); % sum of off-diagonal elements
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
    results.norm_valmat=results.norm_errmat.*vals;
    results.norm_val=sum(results.norm_valmat(:));
end

%% sample inputs
if strcmpi(parser.Results.type,'samp')
    
    % compute outcome counts and error
    [~,samp_counts]=samp_value(dist_1,dist_2,reg,'vals',ones(2),'reg_type',reg_type);
    samp_errcount=samp_value(dist_1,dist_2,reg,'vals',1-eye(2),'reg_type',reg_type);
    samp_err=samp_errcount/sum(samp_counts(:));
    results.samp_errmat=samp_counts;
    results.samp_err=samp_err;
    if optimal_case
        results.samp_dprime=-2*norminv(samp_err); % samp error can be used to compute samp d'
    elseif ~isequal(vals,eye(2)) % if outcome values are supplied
        [samp_val,samp_valmat]=samp_value(dist_1,dist_2,reg,varargin{:});
        results.samp_valmat=samp_valmat;
        results.samp_val=samp_val;
    end    
    
    % Decision variables
    if strcmpi(reg_type,'fun')
        dist_1_cell=num2cell(dist_1,1);
        dv_1=reg(dist_1_cell{:});
        dist_2_cell=num2cell(dist_2,1);
        dv_2=reg(dist_2_cell{:});
    elseif strcmpi(reg_type,'quad')
        f=quad2fun(reg);
        dv_1=f(dist_1')';
        dv_2=f(dist_2')';
    end
    %     f=quad2fun(reg);
    
    %bdv_1=bayes_dec_var(dist_1,mu_1',v_1,mu_2',v_2,priors(1),vals);
    %bdv_2=bayes_dec_var(dist_2,mu_1',v_1,mu_2',v_2,priors(1),vals);
    results.samp_dv={dv_1,dv_2};
    
    %% sample-optimal boundary
    if parser.Results.samp_opt && strcmpi(reg_type,'quad') % && isempty(parser.Results.reg)% if default boundary,
        try
            % find quad boundary that optimizes expected value / accuracy
            x=fminsearch(@(x) -samp_value_flat(x,dist_1,dist_2,vals),[reg.q2(triu(true(size(reg.q2)))); reg.q1(:); reg.q0],optimset('Display','iter','TolX',0,'TolFun',1/(size(dist_1,1)+size(dist_2,1))));
            %             x=fminsearch(@(x) -samp_value_flat(x,dist_1,dist_2,vals),[reg.q2(:); reg.q1(:); reg.q0],optimset('Display','iter','TolFun',1e-6));
            
            q2=zeros(dim);
            q2(triu(true(dim)))=x(1:(dim^2+dim)/2);
            q2=q2+triu(q2,1)';
            
            samp_reg_1.q2=q2;%reshape(x(1:dim^2),[dim dim])';
            samp_reg_1.q1=x(end-dim:end-1);
            %             samp_reg_1.q1=x(dim^2+1:dim^2+dim);
            samp_reg_1.q0=x(end);
            results.samp_opt_bd=samp_reg_1;
            
            % flip boundary sign for b
            samp_reg_2=structfun(@uminus,samp_reg_1,'un',0);
            
            % Decision variables with samp-opt classifier
            f=quad2fun(samp_reg_1);
            dv_samp_1=f(dist_1')';
            dv_samp_2=f(dist_2')';
            results.samp_opt_dv={dv_samp_1,dv_samp_2};
            
            if dim<=3
                % boundary points
                [~,samp_bd_pts_1]=int_norm_along_angles(mu_1,v_1,samp_reg_1,'n_bd_pts',n_samp_bd_pts);
                [~,samp_bd_pts_2]=int_norm_along_angles(mu_2,v_2,samp_reg_2,'n_bd_pts',n_samp_bd_pts);
                samp_bd_pts=uniquetol([samp_bd_pts_1,samp_bd_pts_2]',1e-12,'Byrows',true,'Datascale',1)'; % trim to unique boundary points
                results.samp_opt_bd_pts=samp_bd_pts;
            end
            
%             samp_opt_success=true;
            
        catch %mException
            %warning(getReport(mException,'extended','hyperlinks','on'))
            warning('Cannot optimize quadratic classifier for the sample (usually due to system limits).')
%             samp_opt_success=false;
            
            %             if dim>3
            %                 % find optimal bdv
            %                 samp_opt_bdv=-fminbnd(@(x) -samp_value_flat(1,[0;1;x],bdv_1,bdv_2,vals),min(bdv_2),max(bdv_1),optimset('Display','iter'));
            %                 %x=fminsearch(@(x) -samp_value_flat(1,[0;1;x],lpr_1,lpr_2,vals),0,optimset('TolFun',1e-6,'Display','iter'));
            %                 results.samp_opt_bdv=samp_opt_bdv;
            %
            %                 samp_reg_1.q2=0;
            %                 samp_reg_1.q1=1;
            %                 samp_reg_1.q0=-samp_opt_bdv;
            %                 dist_1=bdv_1;
            %                 dist_2=bdv_2;
            %             end
        end
        
        % compute outcome counts and error with optimized sample boundary
        [~,samp_opt_counts]=samp_value(dist_1,dist_2,samp_reg_1,'vals',ones(2));
        samp_opt_errcount=samp_value(dist_1,dist_2,samp_reg_1,'vals',1-eye(2));
        samp_opt_err=samp_opt_errcount/sum(samp_counts(:));
        results.samp_opt_errmat=samp_opt_counts;
%         results.samp_opt_errs=samp_opt_counts./sum(samp_opt_counts,2);
        results.samp_opt_err=samp_opt_err;
        
        if optimal_case
            results.samp_opt_dprime=-2*norminv(samp_opt_err); % samp error can be used to compute samp d'
        elseif ~isequal(vals,eye(2)) % if outcome values are supplied
            [samp_opt_val,samp_opt_valmat]=samp_value(dist_1,dist_2,samp_reg_1,'vals',vals);
            results.samp_opt_valmat=samp_opt_valmat;
            results.samp_opt_val=samp_opt_val;
        end
    end
    
    
end

%% Plot samples
if bPlot
    hold on
    
    if dim>3
        if isequal(vals,eye(2))
            xlabel('Bayes decision variable: $\ln \frac{p(N_1 | x)}{p(N_2 | x)}$','interpreter','latex','fontsize',15)
        else
            xlabel('Bayes decision variable: $\ln \frac{ \langle v(N_1 | x) \rangle }{ \langle v(N_2 | x) \rangle}$','interpreter','latex','fontsize',15)
        end
    end
    
    if strcmpi(parser.Results.type,'norm')
        title(sprintf("error = %g",norm_err)) % plot title
    elseif strcmpi(parser.Results.type,'samp')
        
        if dim <=3
            plot_sample(dist_1,priors(1),colors(1,:));
            plot_sample(dist_2,priors(2),colors(2,:));
        elseif dim>3
            % plot sample log posterior ratios
            plot_sample(dv_1,priors(1),colors(1,:))
            plot_sample(dv_2,priors(2),colors(2,:))
        end
        if ~exist('samp_opt_err','var') % if custom boundary
            title(sprintf("error = %g / %g",[norm_err,samp_err])) % plot title
            % don't plot sample boundary
        else
            title(sprintf("error = %g / %g / %g",[norm_err,samp_err,samp_opt_err])) % plot title
            if dim <=3
                plot_boundary(samp_reg_1,dim,'reg_type','quad','plot_type','line','line_color',[0 .7 0]);
                %plot_boundary(samp_reg_1,dim,'reg_type','quad','orig',mu_2,'v',v_2,'plot_type','line','plot_color',[0 .7 0]);
                %                 elseif dim>3 && ~samp_opt_success
                %                     xline(samp_opt_bdv,'color',[0 .7 0]);
%             else
%                 title(sprintf("error = %g / %g",[norm_err,samp_err])) % plot title
            end
        end
        % boundary legends
        legend_marker=line(nan, nan,'color','none','linestyle','none');
        if parser.Results.samp_opt && dim <=3
            legend(legend_marker,'\color[rgb]{0,.7,0}sample-optimized boundary')
            legend box off
        end
    end
    hold off
end