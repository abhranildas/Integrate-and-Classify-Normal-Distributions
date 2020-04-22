function plot_boundary(reg,dim,varargin)
parser = inputParser;
addRequired(parser,'reg',@(x) isnumeric(x)||isstruct(x)|| isa(x,'function_handle'));
addParameter(parser,'reg_type','bd_pts');
addRequired(parser,'dim',@(x) (x==1)||(x==2)||(x==3) );
addParameter(parser,'orig',zeros(dim,1), @isnumeric);
addParameter(parser,'v',eye(dim), @isnumeric);
addParameter(parser,'cheb_reg_span',5);
addParameter(parser,'n_bd_pts',1e3);
addParameter(parser,'plot_color','k');
parse(parser,reg,dim,varargin{:});
reg_type=parser.Results.reg_type;
dim=parser.Results.dim;
orig=parser.Results.orig;
v=parser.Results.v;
cheb_reg_span=parser.Results.cheb_reg_span;
n_bd_pts=parser.Results.n_bd_pts;
plot_color=parser.Results.plot_color;

if strcmpi(reg_type,'cheb')
    if nargin(reg)==1
        fimplicit(reg,[repelem(orig',2)+cheb_reg_span*[-1 1] ylim],'color',plot_color)
    elseif nargin(reg)==2
        %reg2=chebfun2(reg,repelem(orig',2)+cheb_reg_span*[-1 1 -1 1]);
        %contour(reg2, [0 0],'linecolor',plot_color);
        fimplicit(reg,repelem(orig',2)+cheb_reg_span*[-1 1 -1 1],'color',plot_color)
    elseif nargin(reg)==3
        fimplicit3(reg,repelem(orig',2)+cheb_reg_span*[-1 1 -1 1 -1 1],'facecolor',plot_color,'facealpha',0.1,'edgecolor',plot_color,'meshdensity',20)
        %reg3=chebfun3(reg,repelem(orig',2)+cheb_reg_span*[-1 1 -1 1 -1 1]);
        %h=isosurface(reg3,0);
        %set(h,'facecolor',plot_color,'edgecolor',plot_color,'facealpha',0.1,'edgealpha',0.2,'marker','.','markersize',3,'markeredgecolor',plot_color)
        %lighting none
    end
else
    if strcmpi(reg_type,'quad') || strcmpi(reg_type,'ray_scan')
        global bd_pts
        bd_pts=[];
        [~,bd_pts]=prob_bd_angle(orig,v,reg,'reg_type',reg_type,'n_bd_pts',n_bd_pts);
    elseif strcmpi(reg_type,'bd_pts')
        bd_pts=reg;
    end
    
    if dim==1
        for x=bd_pts
            xline(x,'color',plot_color,'linewidth',1);
        end
    elseif dim==2
        plot(bd_pts(1,:),bd_pts(2,:),'.','color',plot_color,'markersize',3)
    elseif dim==3
        grid on
        plot3(bd_pts(1,:),bd_pts(2,:),bd_pts(3,:),'.','color',plot_color,'markersize',3);
    end

end
% Create grid of unit vectors
% if dim==1
%     n=1;
% elseif dim==2
%     dth=pi/n_rays;
%     th=0:dth:pi;
%     n=[cos(th);sin(th)];
% elseif dim==3
%     n=fibonacci_sphere(n_rays);
% end
% 
% [~,x]=reg_ray_scan(n,orig);
% bd_pts_rel=cellfun(@(x,y) x.*y,x,num2cell(n,1),'un',0); % convert to co-ordinates
% bd_pts=horzcat(bd_pts_rel{:})+orig;

