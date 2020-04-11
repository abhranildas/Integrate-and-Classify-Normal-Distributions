function [p_ray,bd_pts_angle]=prob_bd_angle(mu,v,reg,varargin)

global bd_pts

% parse inputs
parser = inputParser;
addRequired(parser,'mu',@isnumeric);
addRequired(parser,'v',@isnumeric);
addRequired(parser,'reg',@(x) isstruct(x)|| isa(x,'function_handle'));
addParameter(parser,'side','normal');
addParameter(parser,'reg_type','quad');
addParameter(parser,'cheb_reg_span',5);
addParameter(parser,'n_bd_pts',1e4);
addParameter(parser,'theta',nan,@isnumeric);
addParameter(parser,'phi',nan,@isnumeric);

parse(parser,mu,v,reg,varargin{:});

reg_type=parser.Results.reg_type;
cheb_reg_span=parser.Results.cheb_reg_span;
n_bd_pts=parser.Results.n_bd_pts;
side=parser.Results.side;
theta=parser.Results.theta;
phi=parser.Results.phi;

dim=length(mu);
if any(strcmp(parser.UsingDefaults,'theta')) && any(strcmp(parser.UsingDefaults,'phi'))
    points=fibonacci_sphere(n_bd_pts);
    [az,el]=cart2sph(points(1,:),points(2,:),points(3,:));
    theta=pi/2-el;
    phi=az;
    if dim==2
        theta=linspace(0,pi,n_bd_pts);
    end
end

if dim==1
    n_z=1;
elseif dim==2
    [x,y]=pol2cart(theta,1);
    n_z=[x;y];
elseif dim==3
    az=phi;
    el=pi/2-theta;
    [x,y,z]=sph2cart(az,el,1);
    n_z=[x(:)';y(:)';z(:)'];
end

% rayscan any region
if strcmp(reg_type,'ray_scan') % ray format region
    % get integral and boundary points using ray method
    reg_rayscan=reg;
else % quad or cheb region
    % get integral and boundary points from the ray-scanned chebfun region using ray method
    reg_rayscan=@(n,orig) ray_scan(reg,n,orig,'reg_type',reg_type,'cheb_reg_span',cheb_reg_span);
end

% unit rays in the original space:
n_x=sqrtm(v)*n_z;
n_x=n_x./vecnorm(n_x);

% initial signs and boundary distances along each direction
[init_sign,x]=reg_rayscan(n_x,mu);

% relative boundary points in original space
rel_bd_pts_angle=cellfun(@(x_ray,n_ray) x_ray.*n_ray, x,num2cell(n_x,1),'un',0);

bd_pts_angle=horzcat(rel_bd_pts_angle{:})+mu;
bd_pts=[bd_pts,bd_pts_angle];

% standard boundary points
std_bd_pts_angle=cellfun(@(bd_pt_rel) sqrtm(v)\bd_pt_rel, rel_bd_pts_angle,'un',0);

% standard boundary distances, sorted
z=cellfun(@(n_ray,bd_pt_std) sort(n_ray'*bd_pt_std), num2cell(n_z,1), std_bd_pts_angle,'un',0);

% total integral
if strcmpi(side,'normal')
    p_ray=cellfun(@(init_sign_ray,z_ray) prob_ray(init_sign_ray,z_ray,dim),num2cell(init_sign),z);
elseif strcmpi(side,'complement')
    p_ray=cellfun(@(init_sign_ray,z_ray) prob_ray(init_sign_ray,z_ray,dim,'complement'),num2cell(init_sign),z);
end

if dim==2
    p_ray=p_ray/pi;
elseif dim==3
    p_ray=reshape(p_ray,size(theta));
    p_ray=p_ray.*sin(theta)/(2*pi);
end
