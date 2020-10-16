function [p_rays,bd_pts_rays]=int_norm_along_rays(mu,v,reg,n_z,varargin)

% parse inputs
parser=inputParser;
parser.KeepUnmatched=true;
addRequired(parser,'mu',@isnumeric);
addRequired(parser,'v',@isnumeric);
addRequired(parser,'reg',@(x) isstruct(x)|| isa(x,'function_handle'));
addParameter(parser,'side','normal');
addParameter(parser,'reg_type','quad');
addParameter(parser,'fun_span',3);
addParameter(parser,'fun_resol',100);
addParameter(parser,'fun_level',0);
addParameter(parser,'n_bd_pts',1e4);

parse(parser,mu,v,reg,varargin{:});

reg_type=parser.Results.reg_type;
side=parser.Results.side;

dim=length(mu);

% rayscan any region
if strcmpi(reg_type,'ray_scan') % ray format region
    reg_s_rayscan=@(n) reg(n,mu,v);
else
    reg_s_rayscan=@(n) ray_scan(reg,n,'mu',mu,'v',v,varargin{:});
end
% initial signs and boundary distances in standardized space
[init_sign,z]=reg_s_rayscan(n_z);

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
