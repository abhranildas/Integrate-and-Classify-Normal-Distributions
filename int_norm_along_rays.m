function [p_rays,bd_pts_rays]=int_norm_along_rays(mu,v,dom,n_z,varargin)

% parse inputs
parser=inputParser;
parser.KeepUnmatched=true;
addRequired(parser,'mu',@isnumeric);
addRequired(parser,'v',@isnumeric);
addRequired(parser,'dom',@(x) isstruct(x)|| isa(x,'function_handle'));
addParameter(parser,'side','normal');
addParameter(parser,'dom_type','quad');
addParameter(parser,'fun_span',5);
addParameter(parser,'fun_resol',100);
addParameter(parser,'n_bd_pts',1e4);

parse(parser,mu,v,dom,varargin{:});

side=parser.Results.side;

dim=length(mu);

% ray-trace the domain
dom_standard_raytrace=@(n) standard_ray_trace(dom,n,'mu',mu,'v',v,varargin{:});

% initial signs and boundary distances in standardized space
[init_sign,z]=dom_standard_raytrace(n_z);

% standard boundary points
std_bd_pts_ray=cellfun(@(z_ray,n_ray) z_ray.*n_ray, z,num2cell(n_z,1),'un',0);

% boundary points
bd_pts_rays=sqrtm(v)*horzcat(std_bd_pts_ray{:})+mu;
global bd_pts
bd_pts=[bd_pts,bd_pts_rays];

% probability on ray
if strcmpi(side,'normal')
    p_rays=cellfun(@(init_sign_ray,z_ray) prob_ray(init_sign_ray,z_ray,dim),num2cell(init_sign),z);
elseif strcmpi(side,'complement')
    p_rays=cellfun(@(init_sign_ray,z_ray) prob_ray(init_sign_ray,z_ray,dim,'complement'),num2cell(init_sign),z);
end
