function [p,pc,bd_pts]=integrate_normal(mu,v,reg,varargin)
% Integrate a normal distribution over a specified region.
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
addRequired(parser,'mu',@isnumeric);
addRequired(parser,'v',@isnumeric);
addRequired(parser,'reg',@(x) isstruct(x)|| isa(x,'function_handle'));
addParameter(parser,'reg_type','quad');
addParameter(parser,'fun_span',3); % compute boundary around the mean within 5*semi-major axis of error ellipse
addParameter(parser,'fun_resol',100);
addParameter(parser,'prior',1,@isnumeric);
addParameter(parser,'AbsTol',1e-10);
addParameter(parser,'RelTol',1e-2);
%addParameter(parser,'n_bd_pts',1e4);
addParameter(parser,'bPlot',2);
colors=colororder;
blue=colors(1,:);
addParameter(parser,'plot_color',blue);

parse(parser,mu,v,reg,varargin{:});
reg=parser.Results.reg;
reg_type=parser.Results.reg_type;
fun_span=parser.Results.fun_span;
fun_resol=parser.Results.fun_resol;
prior=parser.Results.prior;
AbsTol=parser.Results.AbsTol;
RelTol=parser.Results.RelTol;
bPlot=parser.Results.bPlot;
plot_color=parser.Results.plot_color;
%n_bd_pts=parser.Results.n_bd_pts;

dim=length(mu);

if dim <=3
    [p,pc,bd_pts]=int_norm_ray(mu,v,reg,'reg_type',reg_type,'fun_span',fun_span,'fun_resol',fun_resol,'AbsTol',AbsTol,'RelTol',RelTol);
    
    % boundary points
    %[~,bd_pts]=prob_bd_angle(mu,v,reg,'reg_type',reg_type,'func_span',func_span,'n_bd_pts',n_bd_pts);
    
    % plot
    if bPlot
        plot_normal(mu,v,prior,plot_color)
        hold on
        if bPlot==2
            plot_boundary(reg,dim,'reg_type',reg_type,'mu',mu,'v',v,'fill_colors',plot_color);
            %         plot_color_faded=mean([1,1,1; plot_color]);
            %         colormap([1 1 1; plot_color_faded]);
        end
        plot_boundary(bd_pts,dim,'reg_type','bd_pts');
        title(sprintf('$p=%g, \\, \\overline{p}=%g$',[p,pc]),'interpreter','latex')
%         title(sprintf("p = %g, p_c = %g",[p,pc])) % plot title
        hold off;
    end
    
elseif strcmpi(reg_type,'quad')
    % get integral from the generalized chi-squared method
    [p,pc]=int_norm_quad_gx2(mu,v,reg,'AbsTol',AbsTol,'RelTol',RelTol);
    bd_pts=[];
    
    if bPlot
        
        % plot distribution of q(x)
        if nnz(reg.q2) % if the quadratic term exists
            % q(x) ~ generalized chi-squared
            [lambda,m,delta,sigma,c]=gx2_params_norm_quad(mu,v,reg);
            [mu_q,v_q]=gx2stat(lambda,m,delta,sigma,c); % mean and variance of q(x)
            %fplot(@(x) prior*gx2pdf(x,lambda,m,delta,sigma,c), 'color',plot_color,'linewidth',1);
            %hold on;
            x=linspace(mu_q-5*sqrt(v_q),mu_q+5*sqrt(v_q),1e3);
            y=prior*arrayfun(@(x) gx2pdf(x,lambda,m,delta,sigma,c),x);
            area(x,y,'facecolor',plot_color,'facealpha',0.4,'edgecolor',plot_color,'edgealpha',0.5,'linewidth',1)
            %hold off
            
        else % q(x) ~ normal
            mu_q=reg.q1'*mu+reg.q0;
            v_q=reg.q1'*v*reg.q1;
            plot_normal(mu_q,v_q,prior,plot_color);
        end
        
        % plot boundary
        if bPlot==2
            hold on;
            plot_boundary(reg,dim,'fill_colors',plot_color);
        end
        
        title(sprintf('$p=%g, \\, \\overline{p}=%g$',[p,pc]),'interpreter','latex')
        xlabel('$q(${\boldmath$x$}$)$','interpreter','latex');
        hold off;
    end
end

