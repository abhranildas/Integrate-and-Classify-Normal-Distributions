function results=classify_normals(dist_1,dist_2,varargin)
    % CLASSIFY_NORMALS Compute quantities concerning classification accuracy
    % between two uni- or multi-variate normal distributions, or samples.
    %
    % Abhranil Das <abhranil.das@utexas.edu>
    % Center for Perceptual Systems, University of Texas at Austin
    % If you use this code, please cite:
    % <a href="matlab:web('https://arxiv.org/abs/2012.14331')"
    % >A method to integrate and classify normal distributions</a>.
    %
    % Example:
    % mu_1=[4; 5]; v_1=[2 1; 1 1];
    % mu_2=mu_1; v_2=3*[2 -1; -1 1];
    % results=classify_normals([mu_1,v_1],[mu_2,v_2])
    %
    % Required inputs:
    % dist_1        normal 1 parameters [mu_1,v_1] (where mu_1 is a column
    %               vector, and v_1 is the covariance matrix), or sample 1,
    %               where columns are variables and rows are observations.
    % dist_2        similar normal 2 parameters [mu_2,v_2], or sample 2.
    %
    % Optional name-value inputs:
    % prior_1       prior probability of distribution 1. Default is 0.5.
    %               For sample inputs, priors are by default the relative
    %               sample sizes. If prior is specified then, it is taken to be
    %               the prior of the normal fitted to the data, and affects
    %               only the normal-fit-based outputs, not the sample-based
    %               outputs.
    % vals          matrix of outcome values. v(i,j) is the value of
    %               classifying a sample from i as j.
    % input_type    'norm' for normal parameter inputs (default), 'samp' for
    %               sample inputs.
    % dom           custom (non-optimal) classification domain, in one of three
    %               forms:
    %               • struct containing coefficients a2 (matrix), a1 (column
    %                 vector) and a0 (scalar) of a quadratic domain:
    %                 x'*a2*x + a1'*x + a0 > 0
    %               • handle to a ray-trace function, returning the starting sign
    %                 and roots of the domain along any ray
    %               • handle to an implicit function f(x) defining the domain f(x)>0.
    % dom_type      'quad' (default), 'ray_trace' or 'fun' for the above three
    %               types resp.
    % fun_span      trace radius (in Mahalanobis distance) for implicit function
    %               domains. Default=5.
    % fun_resol     resolution of tracing (finding roots) of implicit domain.
    %               Default=100.
    % method        Integration method. 'ray' for ray-trace, or 'gx2'
    %               for generalized chi-square (quad domains only). 'ray'
    %               is default upto 3D, and for non-quad domains beyond 3D.
    %               'gx2' is default for quad domains beyond 3D.
    % samp_opt      true (default) if boundary will be optimized for the
    %               sample, otherwise false.
    % AbsTol        absolute tolerance for the error rate computations. Default=1e-10.
    % RelTol        relative tolerance for the error rate computations. Default=1e-2.
    %               The absolute OR the relative tolerance will be satisfied.
    %               They are not used if the no. of dimensions is >3 and
    %               the domain is not a quadratic. Use mc_samples instead.
    % mc_samples    No. of Monte-Carlo samples of rays to compute error rates.
    %               Used only if the no. of dimensions is >3 and the domain is
    %               not a quadratic. Default=500.
    % n_samp_bd_pts number of sample boundary points to be computed
    % plotmode      'norm_prob' (default): normal probability picture, i.e.
    %               plot of the normals and their classification domains,
    %               'fun_prob': function probability picture, i.e. plot of
    %               the 1d pdfs of the scalar function of the normals that
    %               defines the domain. For >3 dimensions, only fun_prob is
    %               possible. For ray-trace domains, only norm_prob is
    %               possible.
    %               false or 0, for no plot.
    %
    % Output: struct containing
    % norm_bd       struct of coefficients of the optimal quadratic boundary
    %               between the normals supplied, or fitted to the samples
    %               supplied.
    % norm_bd_pts   points on the above boundary computed by the ray-trace
    %               integration method.
    % norm_errmat   error matrix. e(i,j)=prob. of classifying a sample from
    %               normal i as j.
    % norm_err      overall error rate
    % norm_d_b      Bayes-optimal discriminability index between the
    %               distributions. Equal to the separation between two
    %               unit-variance normals that have the same overlap.
    % norm_d_a      Simpson and Fitter's discriminability index. Mahalanobis
    %               distance between the normals, that averages the two
    %               covariance matrices.
    % norm_d_e      Egan and Clarke's discriminability index. Mahalanobis
    %               distance between the normals, that averages the two
    %               square roots of the covariance matrices (i.e. the sd's).
    % norm_valmat   matrix of expected classification outcome values. v(i,j)=
    %               expected value of classifying sample from i as j. Returned
    %               only when custom values are supplied.
    % norm_val      overall expected outcome value. Returned only when custom
    %               values are suppliued.
    % samp_errmat   error matrix of supplied samples classified using norm_bd.
    %               e(i,j)= no. of i samples classified as j.
    % samp_err      overall error rate of classifying the samples, using the
    %               optimal boundary between the normals fitted to the samples.
    % samp_d_b      Bayes-optimal discriminability using samp_err.
    %               Returned only when sample sizes are equal, values are
    %               default, and the classifier is optimal.
    % samp_ov_err   overall error rate of classifying the samples, computed using
    %               the overlap area of the two histograms. This is useful
    %               when the optimal classification boundary is unknown,
    %               and the distributions are quite non-normal.
    % samp_ov_d_b   Bayes-optimal discriminability using samp_ov_err.
    %               Returned only when sample sizes are equal, values are
    %               default, and the classifier is not custom.
    % samp_dv       scalar decision variables that the supplied samples are
    %               mapped to using norm_bd
    % samp_valmat   matrix of total sample classification outcome values using
    %               norm_bd. Returned only when custom values are supplied.
    % samp_val      total value of classifying the samples using norm_bd.
    %               Returned only when custom values are supplied.
    % samp_opt_bd   struct containing coefficients of the boundary optimized
    %               for the samples
    % samp_opt_bd_pts   points on samp_opt_bd computed by the ray-trace
    %                   integration method
    % samp_opt_errmat   error matrix (counts) of supplied samples using
    %                   samp_opt_bd.
    % samp_opt_err  overall error rate of classifying the samples using
    %               samp_opt_bd
    % samp_opt_d_b  Bayes-optimal discriminability based on samp_opt_err.
    %               Returned only when sample sizes are equal, values are
    %               default, and the classifier is optimal.
    % samp_opt_dv   scalar decision variables that the supplied samples are
    %               mapped to using samp_opt_bd
    % samp_opt_valmat   matrix of total sample classification outcome values
    %                   using samp_opt_bd. Returned only when custom values are
    %                   supplied.
    % samp_opt_val      total value of classifying the samples using
    %                   samp_opt_bd. Returned only when custom values are
    %                   supplied.
    %
    % See also:
    % <a href="matlab:open(strcat(fileparts(which('integrate_normal')),filesep,'doc',filesep,'GettingStarted.mlx'))">Interactive demos</a>
    % classify_normals_multi
    % integrate_normal

    %% parse inputs
    parser=inputParser;
    parser.KeepUnmatched=true;
    addRequired(parser,'dist_1',@isnumeric);
    addRequired(parser,'dist_2',@isnumeric);
    addParameter(parser,'prior_1',0.5, @(x) isnumeric(x) && isscalar(x) && (x > 0) && (x < 1));
    addParameter(parser,'vals',eye(2), @(x) isnumeric(x) && ismatrix(x));
    addParameter(parser,'dom',[]);
    addParameter(parser,'dom_type','quad');
    addParameter(parser,'method','ray');
    addParameter(parser,'fun_span',5);
    addParameter(parser,'fun_resol',100);
    addParameter(parser,'input_type','norm', @(s) strcmpi(s,'norm') || strcmpi(s,'samp'));
    addParameter(parser,'samp_opt',true, @islogical);
    addParameter(parser,'n_samp_bd_pts',1e4);
    addParameter(parser,'plotmode','norm_prob',@(s) strcmpi(s,'norm_prob') || strcmpi(s,'fun_prob') || s==false);

    parse(parser,dist_1,dist_2,varargin{:});
    dom=parser.Results.dom;
    dom_type=parser.Results.dom_type;
    method=parser.Results.method;
    vals=parser.Results.vals;
    n_samp_bd_pts=parser.Results.n_samp_bd_pts;
    plotmode=parser.Results.plotmode;

    if strcmpi(parser.Results.input_type,'norm')
        mu_1=dist_1(:,1);
        v_1=dist_1(:,2:end);
        mu_2=dist_2(:,1);
        v_2=dist_2(:,2:end);
    elseif strcmpi(parser.Results.input_type,'samp')
        mu_1=mean(dist_1)';
        v_1=cov(dist_1);
        mu_2=mean(dist_2)';
        v_2=cov(dist_2);
    end

    % make covariances symmetric and positive-definite if needed:
    [~,flag]=chol(v_1);
    if flag, v_1=nearestSPD(v_1); end

    [~,flag]=chol(v_2);
    if flag, v_2=nearestSPD(v_2); end

    % if input samp and prior is not specified,
    if strcmpi(parser.Results.input_type,'samp') && any(strcmpi(parser.UsingDefaults,'prior_1'))
        % set prior according to input sample sizes
        priors(1)=size(dist_1,1)/(size(dist_1,1)+size(dist_2,1));
    else
        priors(1)=parser.Results.prior_1;
    end
    priors(2)=1-priors(1);

    % check if optimal case
    if priors(1)==0.5 && isequal(vals,eye(2)) && isempty(dom)
        optimal_case=true;
    else
        optimal_case=false;
    end

    dim=length(mu_1); % dimension

    colors=colororder;
    if strcmpi('dom_type','ray_trace') && dim>3
        plotmode=false;
    end
    if ~isequal(plotmode,false)
        figure; hold on
        if dim>3
            plotmode='fun_prob';
        elseif strcmpi('dom_type','ray_trace')
            plotmode='norm_prob';
        end
    end


    %% normal inputs

    norm_bd_pts=[];
    if isequal(mu_1,mu_2) && isequal(v_1,v_2)
        % if the dists are identical:
        norm_err_1=0.5; norm_acc_1=0.5;
        norm_err_2=0.5; norm_acc_2=0.5;
    else
        % compute accuracy and boundary for each normal, and combine
        if isempty(dom) % optimal domain
            % compute optimal boundary coefficients
            dom=opt_class_quad([mu_1,v_1],[mu_2,v_2],'vals',vals,'prior_1',priors(1));
            results.norm_bd=dom;
        end

        if strcmpi(dom_type,'ray_trace') % ray-traced region functions
            [norm_acc_1,norm_err_1,norm_bd_pts_1]=integrate_normal(mu_1,v_1,dom,'prior',priors(1),'plot_color',colors(1,:),varargin{:});
            dom_inv=@(n,orig) invert_ray_trace_dom(dom,n,'orig',orig);
            [norm_acc_2,norm_err_2,norm_bd_pts_2]=integrate_normal(mu_2,v_2,dom_inv,'prior',priors(2),'plot_color',colors(2,:),varargin{:});
        else
            if strcmpi(dom_type,'quad') % quadratic domain
                [norm_acc_1,norm_err_1,norm_bd_pts_1]=integrate_normal(mu_1,v_1,dom,'prior',priors(1),'plot_color',colors(1,:),varargin{:});
                [norm_err_2,norm_acc_2,norm_bd_pts_2]=integrate_normal(mu_2,v_2,dom,'prior',priors(2),'plot_color',colors(2,:),varargin{:});
            elseif strcmpi(dom_type,'fun') % function domain
                [norm_acc_1,norm_err_1,norm_bd_pts_1]=integrate_normal(mu_1,v_1,dom,'prior',priors(1),'plot_color',colors(1,:),varargin{:});
                if strcmpi(plotmode,'norm_prob'), hold on, plot_boundary(norm_bd_pts_1,dim,'dom_type','bd_pts'), end
                [norm_err_2,norm_acc_2,norm_bd_pts_2]=integrate_normal(mu_2,v_2,dom,'prior',priors(2),'plot_color',colors(2,:),varargin{:});
                if strcmpi(plotmode,'norm_prob'), hold on, plot_boundary(norm_bd_pts_2,dim,'dom_type','bd_pts'), end
            end
            % plot boundary
            if strcmpi(plotmode,'norm_prob')
                plot_boundary(dom,dim,'dom_type',dom_type);
                if strcmpi(method,'gx2')
                    plot_boundary(dom,dim,'dom_type',dom_type,'plot_type','line');
                end
            elseif strcmpi(plotmode,'fun_prob')
                plot_boundary(@(x) x,1,'dom_type','fun','plot_type','fill');
                plot_boundary(@(x) x,1,'dom_type','fun','plot_type','line');
            end
        end
        norm_bd_pts=uniquetol([norm_bd_pts_1,norm_bd_pts_2]',1e-12,'Byrows',true,'Datascale',1)'; % trim to unique boundary points
    end

    if ~isempty(norm_bd_pts)
        results.norm_bd_pts=norm_bd_pts;
    end

    norm_errmat=[[norm_acc_1, norm_err_1]*priors(1); [norm_err_2, norm_acc_2]*priors(2)];
    results.norm_errmat=norm_errmat;
    norm_err=sum(norm_errmat(~eye(2))); % sum of off-diagonal elements
    results.norm_err=norm_err;

    % d'_b
    if optimal_case
        norm_d_b=-2*norminv(norm_err);
    else
        results_opt=classify_normals([mu_1,v_1],[mu_2,v_2],'plotmode',0);
        norm_d_b=results_opt.norm_d_b;
    end
    results.norm_d_b=norm_d_b;

    % d'_a
    s_rms=sqrtm((v_1+v_2)/2);
    results.norm_d_a=norm(s_rms\(mu_1-mu_2));

    % d'_e
    s_avg=(sqrtm(v_1)+sqrtm(v_2))/2;
    results.norm_d_e=norm(s_avg\(mu_1-mu_2));
    % if d'_b is too large to compute in 1d, supply d'_e instead.
    if dim==1 && isinf(norm_d_b)
        results.norm_d_b=results.norm_d_e;
    end

    if ~isequal(vals,eye(2)) % if outcome values are supplied
        results.norm_valmat=results.norm_errmat.*vals;
        results.norm_val=sum(results.norm_valmat(:));
    end

    %% sample inputs
    if strcmpi(parser.Results.input_type,'samp')

        % compute outcome counts and error
        [~,samp_countmat]=samp_value(dist_1,dist_2,dom,'vals',ones(2),'dom_type',dom_type);
        samp_errcount=samp_value(dist_1,dist_2,dom,'vals',1-eye(2),'dom_type',dom_type);
        samp_err=samp_errcount/sum(samp_countmat(:));
        results.samp_errmat=samp_countmat;
        results.samp_err=samp_err;

        % histogram overlap error:
        samp_ov_err=histogram_overlap_err(dist_1,dist_2);
        results.samp_ov_err=samp_ov_err;
        
        if optimal_case
            results.samp_d_b=-2*norminv(samp_err); % samp error can be used to compute samp d'
            results.samp_ov_d_b=-2*norminv(samp_ov_err); % samp error can be used to compute samp d'
        elseif ~isequal(vals,eye(2)) % if outcome values are supplied
            [samp_val,samp_valmat]=samp_value(dist_1,dist_2,dom,varargin{:});
            results.samp_valmat=samp_valmat;
            results.samp_val=samp_val;
        end

        % Decision variables
        if strcmpi(dom_type,'fun')
            dist_1_cell=num2cell(dist_1,1);
            dv_1=dom(dist_1_cell{:});
            dist_2_cell=num2cell(dist_2,1);
            dv_2=dom(dist_2_cell{:});
            results.samp_dv={dv_1,dv_2};
        elseif strcmpi(dom_type,'quad')
            f=quad2fun(dom);
            dv_1=f(dist_1')';
            dv_2=f(dist_2')';
            results.samp_dv={dv_1,dv_2};
        end

        %% sample-optimal boundary
        if parser.Results.samp_opt && strcmpi(dom_type,'quad') % if default boundary,
            try
                % find quad boundary that optimizes expected value / accuracy
                x=fminsearch(@(x) -samp_value_flat(x,dist_1,dist_2,vals),[dom.q2(triu(true(size(dom.q2)))); dom.q1(:); dom.q0],optimset('Display','iter','TolX',0,'TolFun',1/(size(dist_1,1)+size(dist_2,1))));

                q2=zeros(dim);
                q2(triu(true(dim)))=x(1:(dim^2+dim)/2);
                q2=q2+triu(q2,1)';

                samp_dom_1.q2=q2;
                samp_dom_1.q1=x(end-dim:end-1);
                samp_dom_1.q0=x(end);
                results.samp_opt_bd=samp_dom_1;

                % flip boundary sign for b
                samp_dom_2=structfun(@uminus,samp_dom_1,'un',0);

                % Decision variables with samp-opt classifier
                f=quad2fun(samp_dom_1);
                dv_samp_1=f(dist_1')';
                dv_samp_2=f(dist_2')';
                results.samp_opt_dv={dv_samp_1,dv_samp_2};

                if dim<=3
                    % boundary points
                    [~,samp_bd_pts_1]=int_norm_along_angles(mu_1,v_1,samp_dom_1,'n_bd_pts',n_samp_bd_pts);
                    [~,samp_bd_pts_2]=int_norm_along_angles(mu_2,v_2,samp_dom_2,'n_bd_pts',n_samp_bd_pts);
                    samp_bd_pts=uniquetol([samp_bd_pts_1,samp_bd_pts_2]',1e-12,'Byrows',true,'Datascale',1)'; % trim to unique boundary points
                    results.samp_opt_bd_pts=samp_bd_pts;
                end

            catch %mException
                warning('Cannot optimize quadratic classifier for the sample (usually due to system limits).')
            end

            % compute outcome counts and error with optimized sample boundary
            [~,samp_opt_counts]=samp_value(dist_1,dist_2,samp_dom_1,'vals',ones(2));
            samp_opt_errcount=samp_value(dist_1,dist_2,samp_dom_1,'vals',1-eye(2));
            samp_opt_err=samp_opt_errcount/sum(samp_countmat(:));
            results.samp_opt_errmat=samp_opt_counts;
            results.samp_opt_err=samp_opt_err;

            if optimal_case
                results.samp_opt_d_b=-2*norminv(samp_opt_err); % samp error can be used to compute samp d'
            elseif ~isequal(vals,eye(2)) % if outcome values are supplied
                [samp_opt_val,samp_opt_valmat]=samp_value(dist_1,dist_2,samp_dom_1,'vals',vals);
                results.samp_opt_valmat=samp_opt_valmat;
                results.samp_opt_val=samp_opt_val;
            end
        end


    end

    %% Plot samples
    if ~isequal(plotmode,false)
        if strcmpi(plotmode,'fun_prob')
            if isequal(vals,eye(2))
                xlabel('Bayes decision variable: $\ln \frac{p(N_1 | x)}{p(N_2 | x)}$','interpreter','latex','fontsize',15)
            else
                xlabel('Bayes decision variable: $\ln \frac{ \langle v(N_1 | x) \rangle }{ \langle v(N_2 | x) \rangle}$','interpreter','latex','fontsize',15)
            end
        end

        if strcmpi(parser.Results.input_type,'norm')
            title(sprintf("$p_e = %g$",norm_err),'interpreter','latex') % plot title
        elseif strcmpi(parser.Results.input_type,'samp')

            if strcmpi(plotmode,'norm_prob')
                plot_sample(dist_1,priors(1),colors(1,:));
                plot_sample(dist_2,priors(2),colors(2,:));
            elseif strcmpi(plotmode,'fun_prob')
                % plot sample log posterior ratios
                plot_sample(dv_1,priors(1),colors(1,:))
                plot_sample(dv_2,priors(2),colors(2,:))
            end
            if ~exist('samp_opt_err','var') % if custom boundary
                title(sprintf("$p_e = %g / %g$",[norm_err,samp_err]),'interpreter','latex') % plot title
                % don't plot sample boundary
            else
                title(sprintf("$p_e = %g / %g / %g$",[norm_err,samp_err,samp_opt_err]),'interpreter','latex') % plot title
                if strcmpi(plotmode,'norm_prob')
                    plot_boundary(samp_dom_1,dim,'dom_type','quad','plot_type','line','line_color',[0 .7 0]);
                end
            end
            % boundary legends
            legend_marker=line(nan, nan,'color','none','linestyle','none');
            if parser.Results.samp_opt && strcmpi(plotmode,'norm_prob')
                legend(legend_marker,'\color[rgb]{0,.7,0}sample-optimized boundary')
                legend box off
            end
        end
        hold off
    end