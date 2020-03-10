function [p,pc,bd_pts]=integrate_normal(mu,v,varargin)
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
addRequired(parser,'mu',@(x) isnumeric(x));
addRequired(parser,'v',@(x) isnumeric(x));
addParameter(parser,'bd_coeffs',[]);
addParameter(parser,'bd_fn',[]);
addParameter(parser,'n_points',1e4);
addParameter(parser,'p_prior',1, @(x) isnumeric(x));
addParameter(parser,'bPlot',1, @(x) islogical(x));
addParameter(parser,'plot_color','blue');

parse(parser,mu,v,varargin{:});

n_points=parser.Results.n_points;
bd_fn=parser.Results.bd_fn;
dim=length(mu);

% Cholesky decomposition of distribution a:
C=chol(v,'lower');

if ~isempty(parser.Results.bd_coeffs) % quadratic coefficients
    bd_coeffs=parser.Results.bd_coeffs;
    % standardize  coefficients
    bd_coeffs_std.a2=C'*bd_coeffs.a2*C;
    bd_coeffs_std.a1=C'*(2*bd_coeffs.a2*mu+bd_coeffs.a1);
    bd_coeffs_std.a0=mu'*bd_coeffs.a2*mu+bd_coeffs.a1'*mu+bd_coeffs.a0;
    % get integral from the generalized chi-squared method
    [p,pc]=int_stdnorm_quad_gx2(bd_coeffs_std);
end

bd_pts=[];
if dim <=3
    if ~isempty(parser.Results.bd_coeffs) % quadratic coefficients
        % get boundary points from the grid method
        [~,~,bd_pts]=int_norm_grid(mu,v,@(n) quad_bd(n,mu,v,bd_coeffs),n_points);
    elseif ~isempty(parser.Results.bd_fn) % boundary function
        % get both integral and boundary points from the grid method
        [p,pc,bd_pts]=int_norm_grid(mu,v,bd_fn,n_points);
    end
end

% plot
if parser.Results.bPlot && dim<=3
    plot_normal(mu,v,bd_pts,parser.Results.p_prior,parser.Results.plot_color)
end