function bd_handle=plot_boundary(dom,dim,varargin)
    parser = inputParser;
    parser.KeepUnmatched=true;
    addRequired(parser,'dom',@(x) isnumeric(x)||isstruct(x)|| isa(x,'function_handle'));
    addRequired(parser,'dim',@(x) ismember(x,1:x));
    addOptional(parser,'side','upper',@(x) strcmpi(x,'lower') || strcmpi(x,'upper') );
    addParameter(parser,'dom_type','quad');
    addParameter(parser,'mu',zeros(dim,1), @isnumeric);
    addParameter(parser,'v',eye(dim), @isnumeric);
    addParameter(parser,'fun_level',0);
    addParameter(parser,'n_bd_pts',1e3);
    addParameter(parser,'plot_type','fill');
    addParameter(parser,'line_color',[0 0 0]);
    colors=colororder;
    addParameter(parser,'fill_colors',flipud(.1*colors(1:2,:)+.9*[1 1 1]));
    
    parse(parser,dom,dim,varargin{:});
    side=parser.Results.side;
    dom_type=parser.Results.dom_type;
    mu=parser.Results.mu;
    v=parser.Results.v;
    fun_level=parser.Results.fun_level;
    plot_type=parser.Results.plot_type;
    line_color=parser.Results.line_color;
    fill_colors=parser.Results.fill_colors;
    if isrow(fill_colors)
        fill_colors=[1 1 1; .1*fill_colors+.9*[1 1 1]];
    end
    
%     if dim>3
%         plot_boundary(@(x) x-fun_level,1,'dom_type','fun','plot_type','line');
%         plot_boundary(@(x) x-fun_level,1,'dom_type','fun','plot_type','fill','fill_colors',fill_colors);
%     else
        if strcmpi(dom_type,'fun') || strcmpi(dom_type,'quad')
            
            if strcmpi(dom_type,'fun')
                if dim==1
                    f=@(x) dom(x)-fun_level;
                elseif dim==2
                    f=@(x,y) dom(x,y)-fun_level;
                elseif dim==3
                    f=@(x,y,z) dom(x,y,z)-fun_level;
                end
            elseif strcmpi(dom_type,'quad')
                f=quad2fun(dom,fun_level,1);
            end
            
            if dim==1 || dim==2
                xl=xlim; yl=ylim;
                if strcmpi(plot_type,'line')
                    bd_handle=fimplicit(f,'color',line_color,'linewidth',.75);
                elseif strcmpi(plot_type,'fill')
                    fh=fcontour(f,'Fill','on');
                    clim([-eps eps])
                    if strcmpi(side,'upper')
                        colormap(fill_colors)
                    elseif strcmpi(side,'lower')
                        colormap(flipud(fill_colors))
                    end
                    uistack(fh,'bottom')
                end
                xlim(xl); ylim(yl);
            elseif dim==3
                fimplicit3(f,'facecolor',line_color,'facealpha',0.1,'edgecolor','none','meshdensity',20)
            end
            
        elseif strcmpi(dom_type,'ray_trace') || strcmpi(dom_type,'bd_pts')
            
            if strcmpi(dom_type,'ray_trace')
                global bd_pts
                bd_pts=[];
                [~,bd_pts]=norm_prob_across_angles(mu,v,dom,varargin{:},'bd_pts',true);
            elseif strcmpi(dom_type,'bd_pts')
                bd_pts=dom;
            end
            
            if dim==1
                for x=bd_pts
                    xline(x,'color',line_color,'linewidth',1);
                end
            elseif dim==2
                if ~isempty(bd_pts)
                    plot(bd_pts(1,:),bd_pts(2,:),'.','color',line_color,'markersize',3)
                end
            elseif dim==3
                grid on
                if ~isempty(bd_pts)
                    plot3(bd_pts(1,:),bd_pts(2,:),bd_pts(3,:),'.','color',line_color,'markersize',3);
                end
            end
            
        end
%     end