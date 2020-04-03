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
addParameter(parser,'n_rays',1e4);
addParameter(parser,'prior',1,@isnumeric);
addParameter(parser,'estimate',[],@(x) (isnumeric(x)&&isscalar(x))||strcmpi(x,'tail'));
addParameter(parser,'bPlot',1,@islogical);
colors=colororder;
color=colors(1,:);
addParameter(parser,'plot_color',color);

parse(parser,mu,v,reg,varargin{:});

n_rays=parser.Results.n_rays;
dim=length(mu);

if strcmp(parser.Results.reg_type,'quad') % if quadratic coefficients supplied
    % standardize  coefficients
    reg_coeffs_std.a2=sqrtm(v)*reg.a2*sqrtm(v);
    reg_coeffs_std.a1=sqrtm(v)*(2*reg.a2*mu+reg.a1);
    reg_coeffs_std.a0=mu'*reg.a2*mu+reg.a1'*mu+reg.a0;
    % get integral from the generalized chi-squared method
    if isempty(parser.Results.estimate)
        [p,pc]=int_stdnorm_quad_gx2(reg_coeffs_std);
    else
        [p,pc]=int_stdnorm_quad_gx2(reg_coeffs_std,parser.Results.estimate);
    end
end

bd_pts=[];
if dim <=3
    if strcmp(parser.Results.reg_type,'quad') % if quadratic coefficients supplied
        % get boundary points from the grid method
        [~,~,bd_pts]=int_norm_grid(mu,v,@(n) ray_scan(reg,'quad',n,mu),n_rays);
    elseif strcmp(parser.Results.reg_type,'ray_scan') % ray format region
        % get both integral and boundary points using ray method
        [p,pc,bd_pts]=int_norm_grid(mu,v,reg,n_rays);
    elseif strcmp(parser.Results.reg_type,'cheb') % chebfun region
        % get both integral and boundary points from the ray-scanned chebfun region using ray method
        [p,pc,bd_pts]=int_norm_grid(mu,v,@(n) ray_scan(reg,n,mu),n_rays);
    end
end

% plot
if parser.Results.bPlot && dim<=3
    plot_normal(mu,v,bd_pts,parser.Results.prior,parser.Results.plot_color)
    title(sprintf("p = %g, p_c = %g",[p,pc])) % plot title
end