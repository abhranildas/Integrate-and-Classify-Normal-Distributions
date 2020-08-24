function plot_boundary(reg,dim,varargin)
parser = inputParser;
addRequired(parser,'reg',@(x) isnumeric(x)||isstruct(x)|| isa(x,'function_handle'));
addParameter(parser,'reg_type','bd_pts');
addRequired(parser,'dim',@(x) ismember(x,1:x));
addParameter(parser,'mu',zeros(dim,1), @isnumeric);
addParameter(parser,'v',eye(dim), @isnumeric);
addParameter(parser,'n_bd_pts',1e3);
addParameter(parser,'plot_type','fill');
addParameter(parser,'line_color',[0 0 0]);
colors=colororder;
addParameter(parser,'fill_colors',flipud(.1*colors(1:2,:)+.9*[1 1 1]));
parse(parser,reg,dim,varargin{:});
reg_type=parser.Results.reg_type;
dim=parser.Results.dim;
mu=parser.Results.mu;
v=parser.Results.v;
n_bd_pts=parser.Results.n_bd_pts;
plot_type=parser.Results.plot_type;
line_color=parser.Results.line_color;
fill_colors=parser.Results.fill_colors;
if isrow(fill_colors)
    fill_colors=[1 1 1; .1*fill_colors+.9*[1 1 1]];
end

if dim>3
    plot_boundary(@(x) x,1,'reg_type','fun','plot_type','line');
    plot_boundary(@(x) x,1,'reg_type','fun','plot_type','fill','fill_colors',fill_colors);
    %xline(0,'color',line_color,'linewidth',1);
else
    if strcmpi(reg_type,'fun') || strcmpi(reg_type,'quad')
        
        if strcmpi(reg_type,'fun')
            f=reg;
        elseif strcmpi(reg_type,'quad')
            f=quad2fun(reg,1);
        end
        
        if dim==1 || dim==2
            xl=xlim; yl= ylim;
            if strcmpi(plot_type,'line')
                fimplicit(f,'color',line_color);
            elseif strcmpi(plot_type,'fill')
                fh=fcontour(f,'Fill','on');
                caxis([-eps eps])
                colormap(fill_colors)
                uistack(fh,'bottom')
            end
            xlim(xl); ylim (yl);
        elseif dim==3
            fimplicit3(f,'facecolor',line_color,'facealpha',0.1,'edgecolor','none','meshdensity',20)
        end
        
    elseif strcmpi(reg_type,'ray_scan') || strcmpi(reg_type,'bd_pts')
        
        if strcmpi(reg_type,'ray_scan')
            global bd_pts
            bd_pts=[];
            [~,bd_pts]=prob_bd_angle(mu,v,reg,'reg_type',reg_type,'n_bd_pts',n_bd_pts);
        elseif strcmpi(reg_type,'bd_pts')
            bd_pts=reg;
        end
        
        if dim==1
            for x=bd_pts
                xline(x,'color',line_color,'linewidth',1);
            end
        elseif dim==2
            plot(bd_pts(1,:),bd_pts(2,:),'.','color',line_color,'markersize',3)
        elseif dim==3
            grid on
            plot3(bd_pts(1,:),bd_pts(2,:),bd_pts(3,:),'.','color',line_color,'markersize',3);
        end
        
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

