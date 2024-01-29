function [init_sign,z,samp_correct]=standard_ray_trace(dom,n,varargin)

    % returns a domain in ray-trace format in the standardized space
    % parse inputs
    dim=size(n,1);
    parser=inputParser;
    parser.KeepUnmatched=true;
    addRequired(parser,'dom',@(x) isstruct(x) || isa(x,'function_handle') || ismatrix(x));
    addRequired(parser,'n',@(x) isnumeric(x) || isa(x,'sym'));
    addParameter(parser,'mu',zeros(dim,1),@isnumeric);
    addParameter(parser,'v',eye(dim),@isnumeric);
    addParameter(parser,'dom_type','quad');
    addParameter(parser,'fun_span',5);
    addParameter(parser,'fun_resol',100);
    addParameter(parser,'fun_level',0);

    parse(parser,dom,n,varargin{:});
    mu=parser.Results.mu;
    v=parser.Results.v;
    dom_type=parser.Results.dom_type;
    fun_span=parser.Results.fun_span;
    fun_resol=parser.Results.fun_resol;
    fun_level=parser.Results.fun_level;

    if strcmpi(dom_type,'ray_trace') || strcmpi(dom_type,'rect')
        n_x=sqrtm(v)*n; % transform rays to original space
        if strcmpi(dom_type,'ray_trace')
            [init_sign,x]=dom(n_x,mu); % trace in the original space from mu
        elseif strcmpi(dom_type,'rect')
            [init_sign,x]=rectangle_ray_trace(dom,n_x,mu); % trace in the original space from mu
        end
        z=cellfun(@(a,b) a/b, x,num2cell(vecnorm(n_x)),'un',0); % scale back to standard space

    elseif strcmpi(dom_type,'quad')

        if nargout==3
            [~,~,samp_correct]=samp_value(n',n',dom);
            init_sign=[];
            z=[];
        else
            if isnumeric(n)
                n=n./vecnorm(n); % normalize direction vectors
            end

            % standardized boundary coefficients
            quad_s=standard_quad(dom,mu,v);

            q2=dot(n,quad_s.q2*n);
            q1=quad_s.q1'*n;
            q0=quad_s.q0;

            % sign of the quadratic at -inf:
            % if isnumeric(n)
                % init_sign=sign(q2); % square term sets the sign
                % init_sign(~init_sign)=-sign(q1(~init_sign)); % linear term sets the sign for leftovers
                % init_sign(~init_sign)=sign(q0);% constant term sets the sign for the leftovers
            % elseif isa(n,'sym')
                % init_sign=piecewise(q2~=0,sign(q2),(q2==0)&&(q1~=0),-sign(q1),(q2^2+q1^2==0),sign(q0));
                init_sign=sign(4*sign(q2)-2*sign(q1)+sign(q0));
            % end

            z=arrayfun(@roots_ray_quad,q2,q1,q0*ones(1,size(n,2)),'un',0); % this allows function to calculate on multiple directions at once

        end

    elseif strcmpi(dom_type,'fun')
        if nargout==3
            [~,~,samp_correct]=samp_value(n',n',dom,'dom_type','fun');
            init_sign=[];
            z=[];
        else
            [init_sign,z]=cellfun(@(n_ray) ray_trace_fun(@(r)...
                standard_ray_fun(dom,mu,v,n_ray,r,fun_level),...
                fun_span*[-1 1],fun_resol),num2cell(n,1),'un',0);
            init_sign=cell2mat(init_sign);
        end
    end

end