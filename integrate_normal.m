function [p,p_err,bd_pts]=integrate_normal(mu,v,dom,varargin)
    % INTEGRATE_NORMAL Integrate a (multi)normal in any domain.
    %
    % Abhranil Das <abhranil.das@utexas.edu>
    % Center for Perceptual Systems, University of Texas at Austin
    % If you use this code, please cite:
    % <a href="matlab:web('https://arxiv.org/abs/2012.14331')"
    % >A method to integrate and classify normal distributions</a>.
    %
    % Example:
    % mu=[-1; -1]; v=[1 0.5; 0.5 2];
    % fun=@(x,y) x.*sin(y) - y.*cos(x) -1;
    % [p,~,bd_pts]=integrate_normal(mu,v,fun,'dom_type','fun')
    % [pc,~,bd_pts]=integrate_normal(mu,v,fun,'upper','dom_type','fun')
    %
    % Required inputs:
    % mu            normal mean as column vector
    % v             normal variance-covariance matrix
    % dom           integration domain, in one of four forms:
    %               • struct containing coefficients a2 (matrix), a1 (column
    %                 vector) and a0 (scalar) of a quadratic domain:
    %                 x'*a2*x + a1'*x + a0 > 0
    %               • 2-row matrix, where the first and second row are
    %                 lower and upper limits of a (hyper-)rectangle
    %               • handle to a ray-trace function, returning the starting sign
    %                 and roots of the domain along any ray
    %               • handle to an implicit function f(x) defining the domain f(x)>0.
    %
    % Optional positional input:
    % side          'upper' (default) to compute normal prob. within domain f(x)>0.
    %               'lower' to compute in the complement domain f(x)<0.
    % 
    % Optional name-value inputs:
    % dom_type      'quad' (default), 'rect', 'ray_trace' or 'fun' for the
    %               above four types resp.
    % method        Integration method. 'ray' (default) for ray-trace, or 'gx2'
    %               for generalized chi-square (quad domains only).
    % force_mc      set to true to force the ray method to use Monte-Carlo 
    %               integration instead of grid integration, even for dimensions <=3.
    % n_rays        No. of Monte-Carlo samples of rays. Default=500.
    % fun_span      trace radius (in Mahalanobis distance) for implicit function
    %               domains. Default=5.
    % fun_resol     resolution of tracing (finding roots) of implicit domain.
    %               Default=100.
    % fun_level     level c for defining domain as f(x)>c. Default=0.
    % bd_pts        default=false. true for also returning and plotting boundary points
    %               computed when using ray method.
    % prior         prior probability. Only used for scaling plots.
    %               Default=1.
    % AbsTol        absolute tolerance for the integral. Default=1e-10.
    % RelTol        relative tolerance for the integral. Default=1e-2.
    %               The absolute OR the relative tolerance will be satisfied.
    %               They are not used if the no. of dimensions is >3 and
    %               the domain is not a quadratic. Use n_rays instead.
    % vpa           false (default) to do ray method or Imhof's method integrals numerically,
    %               true to do them symbolically with variable precision.
    % plotmode      'norm_prob' (default): normal probability picture, i.e.
    %               plot of the normal and the domain,
    %               'fun_prob': function probability picture, i.e. plot of
    %               the 1d pdf of the scalar function of the normal that
    %               defines the domain. For >3 dimensions, only fun_prob is
    %               possible. For ray-trace domains, only norm_prob is
    %               possible.
    %               false or 0, for no plot.
    % plot_color    2-row array of plot colors of the normal and the
    %               domain respectively. Single row to skip coloring the
    %               domain.
    %
    % Outputs:
    % p             integrated probability
    % pc            complement of the probability (more accurate when it is
    %               small)
    % bd_pts        points on the domain boundary computed by the ray-trace
    %               integration method. Returned only when bd_pts=true.
    %
    % See also:
    % <a href="matlab:open(strcat(fileparts(which('integrate_normal')),filesep,'doc',filesep,'GettingStarted.mlx'))">Interactive demos</a>
    % norm_fun_cdf
    % classify_normals

    % parse inputs
    parser=inputParser;
    parser.KeepUnmatched=true;
    addRequired(parser,'mu',@isnumeric);
    addRequired(parser,'v',@isnumeric);
    addRequired(parser,'dom',@(x) isstruct(x) || isa(x,'function_handle') || ismatrix(x));
    addOptional(parser,'side','upper',@(x) strcmpi(x,'lower') || strcmpi(x,'upper') );
    addParameter(parser,'dom_type','quad');
    addParameter(parser,'method','ray');
    addParameter(parser,'bd_pts',false);
    addParameter(parser,'fun_span',5);
    addParameter(parser,'fun_resol',100);
    addParameter(parser,'prior',1,@isnumeric);
    addParameter(parser,'AbsTol',1e-10);
    addParameter(parser,'RelTol',1e-2);
    addParameter(parser,'plotmode','norm_prob');
    colors=colororder;
    addParameter(parser,'plot_color',[colors(1,:);colors(1,:)]);

    parse(parser,mu,v,dom,varargin{:});
    dom_type=parser.Results.dom_type;
    method=parser.Results.method;
    prior=parser.Results.prior;
    plotmode=parser.Results.plotmode;
    plot_color=parser.Results.plot_color;
    bd_pts_flag=parser.Results.bd_pts;


    dim=length(mu);

    % if it's a ray-trace domain, return boundary points
    if strcmpi(dom_type,'ray_trace') && ismember('bd_pts',parser.UsingDefaults) && ~bd_pts_flag
        bd_pts_flag=true;
    end

    if any(strcmpi(parser.UsingDefaults,'method')) && dim>3 && strcmpi(dom_type,'quad')
        method='gx2';
    end

    if strcmpi(method,'ray')
        [p,p_err,bd_pts]=int_norm_ray(mu,v,dom,varargin{:},'bd_pts',bd_pts_flag);
    elseif strcmpi(method,'gx2')
        p=int_norm_quad_gx2(mu,v,dom,varargin{:});
        p_err=[];
        bd_pts=[];
        % if gx2 returns 0, use ray that can use vpa
        if p==0 && any(strcmpi(parser.UsingDefaults,'method'))
            [p,p_err,bd_pts]=int_norm_ray(mu,v,dom,varargin{:});
        end
    end

    % plotting
    if strcmpi(dom_type,'ray_trace') && dim>3
        plotmode=false;
    end
    if ~isequal(plotmode,false)
        holdon=ishold;
        if dim>3
            plotmode='fun_prob';
        elseif strcmpi(dom_type,'ray_trace')
            plotmode='norm_prob';
        end
        if strcmpi(plotmode,'norm_prob')
            plot_normal(mu,v,prior,plot_color(1,:))
            hold on
            if size(plot_color,1)==2
                plot_boundary(dom,dim,varargin{:},'mu',mu,'v',v,'fill_colors',plot_color(2,:));                
            end            
            if strcmpi(method,'ray') && bd_pts_flag
                plot_boundary(bd_pts,dim,'dom_type','bd_pts');
            else
                % if classify_normals called this function, don't draw
                % boundary
                callStack = dbstack;
                if ~(numel(callStack)>1 && strcmpi(callStack(2).name,'classify_normals'))
                    plot_boundary(dom,dim,varargin{:},'mu',mu,'v',v,'plot_type','line');
                end
            end
        elseif strcmpi(plotmode,'fun_prob')
            ylim auto
            plot_norm_fun(mu,v,dom,prior,plot_color(1,:),varargin{:})
            if size(plot_color,1)==2
                hold on
                plot_boundary(@(x) x,1,'dom_type','fun','plot_type','fill','fill_colors',plot_color(2,:));
                plot_boundary(@(x) x,1,'dom_type','fun','plot_type','line');
            end
            % set ylim to start from 0
            yl=ylim;
            ylim([0 yl(2)])
        end
        title(sprintf('$p=%g$',p),'interpreter','latex')
        if ~holdon
            hold off
        end
    end
