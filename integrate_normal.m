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
addParameter(parser,'cheb_reg_span',3); % compute boundary around the mean within 5*semi-major axis of error ellipse
addParameter(parser,'func_crossings',100);
addParameter(parser,'prior',1,@isnumeric);
addParameter(parser,'AbsTol',1e-10);
addParameter(parser,'RelTol',1e-2);
%addParameter(parser,'n_bd_pts',1e4);
addParameter(parser,'bPlot',1,@islogical);
colors=colororder;
color=colors(1,:);
addParameter(parser,'plot_color',color);

parse(parser,mu,v,reg,varargin{:});
reg=parser.Results.reg;
reg_type=parser.Results.reg_type;
cheb_reg_span=parser.Results.cheb_reg_span;
func_crossings=parser.Results.func_crossings;
prior=parser.Results.prior;
AbsTol=parser.Results.AbsTol;
RelTol=parser.Results.RelTol;
bPlot=parser.Results.bPlot;
plot_color=parser.Results.plot_color;
%n_bd_pts=parser.Results.n_bd_pts;

dim=length(mu);

C=sqrtm(v);
[~,D] = eig(C);

if dim <=3
    [p,pc,bd_pts]=int_norm_grid(mu,v,reg,'reg_type',reg_type,'cheb_reg_span',cheb_reg_span*max(diag(D)),'func_crossings',func_crossings,'AbsTol',AbsTol,'RelTol',RelTol);
    
    % boundary points
    %[~,bd_pts]=prob_bd_angle(mu,v,reg,'reg_type',reg_type,'cheb_reg_span',cheb_reg_span,'n_bd_pts',n_bd_pts);
    
    % plot
    if bPlot
        plot_normal(mu,v,prior,plot_color)
        hold on;
        %         if strcmpi(reg_type,'cheb')
        % plot cheb region in grey
        %             plot_boundary(reg,dim,'reg_type',reg_type,'orig',mu,'v',v,'cheb_reg_span',cheb_reg_span,'plot_color',.5*[1 1 1]);
        % plot only the cheb region used for integration in black
        %             plot_boundary(bd_pts,dim,'reg_type','bd_pts');
        %         else
        if dim==1
            plot_boundary(bd_pts,dim,'reg_type','bd_pts');
        else
            plot_boundary(reg,dim,'reg_type',reg_type,'orig',mu,'v',v,'cheb_reg_span',cheb_reg_span*max(diag(D)));
        end
        title(sprintf("p = %g, p_c = %g",[p,pc])) % plot title
        hold off;
    end
    
elseif strcmpi(reg_type,'quad')
    % get integral from the generalized chi-squared method
    [p,pc]=int_norm_quad_gx2(mu,v,reg,'AbsTol',AbsTol,'RelTol',RelTol);
    bd_pts=[];
end

