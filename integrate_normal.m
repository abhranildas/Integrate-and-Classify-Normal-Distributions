function [p,pc,bd_pts]=integrate_normal(mu,v,dom,varargin)
% Integrate a normal distribution over a specified domain.
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
parser=inputParser;
parser.KeepUnmatched=true;
addRequired(parser,'mu',@isnumeric);
addRequired(parser,'v',@isnumeric);
addRequired(parser,'dom',@(x) isstruct(x)|| isa(x,'function_handle'));
addParameter(parser,'dom_type','quad');
addParameter(parser,'method','ray');
addParameter(parser,'fun_span',3);
addParameter(parser,'fun_resol',100);
addParameter(parser,'fun_level',0);
addParameter(parser,'prior',1,@isnumeric);
addParameter(parser,'AbsTol',1e-10);
addParameter(parser,'RelTol',1e-2);
addParameter(parser,'plotmode',2);
colors=colororder;
addParameter(parser,'plot_color',colors(1,:));

parse(parser,mu,v,dom,varargin{:});
dom_type=parser.Results.dom_type;
method=parser.Results.method;
prior=parser.Results.prior;
plotmode=parser.Results.plotmode;
plot_color=parser.Results.plot_color;

dim=length(mu);

if any(strcmpi(parser.UsingDefaults,'method')) && dim>3 && strcmpi(dom_type,'quad')
    method='gx2';
end

if strcmpi(method,'ray')
    [p,pc,bd_pts]=int_norm_ray(mu,v,dom,varargin{:});
elseif strcmpi(method,'gx2')
    % get integral from the generalized chi-squared method
    [p,pc]=int_norm_quad_gx2(mu,v,dom,varargin{:});
    bd_pts=[];
end

% plot
if plotmode
    if dim<=3
        plot_normal(mu,v,prior,plot_color)
        hold on
        if plotmode==2
            plot_boundary(dom,dim,'mu',mu,'v',v,'fill_colors',plot_color,varargin{:});
        end
        if strcmpi(method,'ray')
            plot_boundary(bd_pts,dim,'dom_type','bd_pts');
        end
    else
        if strcmpi(dom_type,'quad')
            % plot distribution of q(x)
            if nnz(dom.q2) % if the quadratic term exists
                % q(x) ~ generalized chi-squared
                [lambda,m,delta,sigma,c]=gx2_params_norm_quad(mu,v,dom);
                [mu_q,v_q]=gx2stat(lambda,m,delta,sigma,c); % mean and variance of q(x)
                x=linspace(mu_q-5*sqrt(v_q),mu_q+5*sqrt(v_q),1e3);
                y=prior*arrayfun(@(x) gx2pdf(x,lambda,m,delta,sigma,c),x);
                area(x,y,'facecolor',plot_color,'facealpha',0.4,'edgecolor',plot_color,'edgealpha',0.5,'linewidth',1)
            else % q(x) ~ normal
                mu_q=dom.q1'*mu+dom.q0;
                v_q=dom.q1'*v*dom.q1;
                plot_normal(mu_q,v_q,prior,plot_color);
            end
            hold on
            % plot boundary
            if plotmode==2
                plot_boundary(dom,dim,'fill_colors',plot_color);
            end
            xlabel('$q(${\boldmath$x$}$)$','interpreter','latex');
        end
    end
    hold off;
    title(sprintf('$p=%g, \\, \\overline{p}=%g$',[p,pc]),'interpreter','latex')
end

end

