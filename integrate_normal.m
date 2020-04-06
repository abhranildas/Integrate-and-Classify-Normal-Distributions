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
addParameter(parser,'prior',1,@isnumeric);
addParameter(parser,'AbsTol',1e-10);
addParameter(parser,'RelTol',1e-6);
addParameter(parser,'bPlot',1,@islogical);
colors=colororder;
color=colors(1,:);
addParameter(parser,'plot_color',color);

parse(parser,mu,v,reg,varargin{:});
AbsTol=parser.Results.AbsTol;
RelTol=parser.Results.RelTol;

dim=length(mu);

% reg rayscan function
if strcmp(parser.Results.reg_type,'quad') % if quadratic coefficients supplied
    % get integral and boundary points from the grid method
    reg_fn=@(n) ray_scan(reg,'quad',n,mu);
elseif strcmp(parser.Results.reg_type,'ray_scan') % ray format region
    % get integral and boundary points using ray method
    reg_fn=reg;
elseif strcmp(parser.Results.reg_type,'cheb') % chebfun region
    % get integral and boundary points from the ray-scanned chebfun region using ray method
    reg_fn=@(n) ray_scan(reg,n,mu);
end

if dim <=3
    [p,pc]=int_norm_grid(mu,v,reg_fn,'AbsTol',AbsTol,'RelTol',RelTol);    
    % boundary points
    n_rays=1e4;
    if dim==1
    [~,bd_pts]=prob_theta(mu,v,reg_fn,nan,nan);
    elseif dim==2
    [~,bd_pts]=prob_theta(mu,v,reg_fn,linspace(0,pi,n_rays),nan);
    elseif dim==3
        points=fibonacci_sphere(n_rays);
        [theta,phi]=cart2sph(points(1,:),points(2,:),points(3,:));
        [~,bd_pts]=prob_theta(mu,v,reg_fn,theta,phi);
    end
elseif strcmp(parser.Results.reg_type,'quad')
    % get integral from the generalized chi-squared method
    [p,pc]=int_norm_quad_gx2(mu,v,reg,'AbsTol',AbsTol,'RelTol',RelTol);
    bd_pts=[];
end

% plot
if parser.Results.bPlot && dim<=3
    plot_normal(mu,v,bd_pts,parser.Results.prior,parser.Results.plot_color)
    title(sprintf("p = %g, p_c = %g",[p,pc])) % plot title
end