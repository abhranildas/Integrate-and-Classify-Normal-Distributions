function plot_norm_fun(mu,v,fun,prior,plot_color,varargin)
    
    if ~exist('prior','var')
        prior=1;
    end
    
    if ~exist('plot_color','var')
        colors=colororder;
        plot_color=colors(1,:);
    end
    
    dim=length(mu);
    
    if isstruct(fun)
        % plot distribution of q(x)
        if nnz(fun.q2) % if the quadratic term exists
            % q(x) ~ generalized chi-squared
            [w,k,lambda,m,s]=gx2_params_norm_quad(mu,v,fun);
            [f,x]=gx2pdf('full',w,k,lambda,m,s);
            area(x,prior*f,'facecolor',plot_color,'facealpha',0.4,'edgecolor',plot_color,'edgealpha',0.5,'linewidth',1)
        else % q(x) ~ normal
            mu_q=fun.q1'*mu+fun.q0;
            v_q=fun.q1'*v*fun.q1;
            plot_normal(mu_q,v_q,prior,plot_color);
        end
        if dim==1
            xlabel('$q(x)$','interpreter','latex');
        else
            xlabel('$q(${\boldmath$x$}$)$','interpreter','latex');
        end
    elseif isa(fun,'function_handle')
        xlims=norm_fun_inv([.01 .99],mu,v,fun,varargin{:});
        x=linspace(xlims(1),xlims(2),100);
        p=norm_fun_cdf(x,mu,v,fun,varargin{:});
        f=diff(p);x=x(1:end-1);
        area(x,prior*f,'facecolor',plot_color,'facealpha',0.4,'edgecolor',plot_color,'edgealpha',0.5,'linewidth',1)
        if dim==1
            xlabel('$f(x)$','interpreter','latex');
        else
            xlabel('$f(${\boldmath$x$}$)$','interpreter','latex');
        end
    end
