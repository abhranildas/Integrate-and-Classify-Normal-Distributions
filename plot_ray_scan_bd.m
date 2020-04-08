function plot_ray_scan_bd(reg_ray_scan,dim,varargin)
parser = inputParser;
addRequired(parser,'reg_ray_scan',@(x) isa(x,'function_handle'));
addRequired(parser,'dim',@(x) (x==1)||(x==2)||(x==3) );
addParameter(parser,'orig',zeros(dim,1), @isnumeric);
addParameter(parser,'n_rays',1e4,@isnumeric);
addParameter(parser,'color','k');
parse(parser,reg_ray_scan,dim,varargin{:});
orig=parser.Results.orig;
n_rays=parser.Results.n_rays;
color=parser.Results.color;

% Create grid of unit vectors
if dim==1
    n=1;
elseif dim==2
    dth=pi/n_rays;
    th=0:dth:pi;
    n=[cos(th);sin(th)];
elseif dim==3
    n=fibonacci_sphere(n_rays);
end

[~,x]=reg_ray_scan(n,orig);
bd_pts_rel=cellfun(@(x,y) x.*y,x,num2cell(n,1),'un',0); % convert to co-ordinates
bd_pts=horzcat(bd_pts_rel{:})+orig;

if dim==1
    for x=bd_pts
        xline(x,'color',color,'linewidth',1);
    end
elseif dim==2
     plot(bd_pts(1,:),bd_pts(2,:),'.','color',color,'markersize',3)
elseif dim==3
    plot3(bd_pts(1,:),bd_pts(2,:),bd_pts(3,:),'.','color',color,'markersize',1);
    grid on
end