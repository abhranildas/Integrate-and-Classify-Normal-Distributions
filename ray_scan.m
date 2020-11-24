function [init_sign,x,samp_correct]=ray_scan(dom,n,varargin)

% parse inputs
dim=size(n,1);
parser=inputParser;
parser.KeepUnmatched=true;
addRequired(parser,'dom',@(x) isstruct(x)|| isa(x,'function_handle'));
addRequired(parser,'n',@isnumeric);
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

% root(s) along a ray through quad domain
    function x_ray=roots_ray_quad(q2_pt,q1_pt)
        x_ray=sort(roots([q2_pt q1_pt q0]))';
        x_ray=x_ray(~imag(x_ray)); % only real roots
        
        % remove any roots that are tangents. Only crossing points
        slope=2*q2_pt*x_ray+q1_pt;
        x_ray=x_ray(slope~=0);
    end

if strcmpi(dom_type,'quad')
    
    if nargout==3
        [~,~,samp_correct]=samp_value(n',n',dom);
        init_sign=[];
        x=[];
    else
        n=n./vecnorm(n); % normalize direction vectors
        
        % standardized boundary coefficients
        quad_s=standard_quad(dom,mu,v);
        
        q2=dot(n,quad_s.q2*n);
        q1=quad_s.q1'*n;
        q0=quad_s.q0;
        
        % sign of the quadratic at -inf:
        init_sign=sign(q2); % square term sets the sign
        init_sign(~init_sign)=-sign(q1(~init_sign)); % linear term sets the sign for leftovers
        init_sign(~init_sign)=sign(q0);% constant term sets the sign for the leftovers
        
        x=arrayfun(@roots_ray_quad,q2,q1,'un',0); % this allows function to calculate on multiple directions at once
    end
    
elseif strcmpi(dom_type,'fun')
    if nargout==3
        [~,~,samp_correct]=samp_value(n',n',dom,'dom_type','fun');
        init_sign=[];
        x=[];
    else
        [init_sign,x]=cellfun(@(n_ray) ray_scan_fun(@(r)...
            standard_ray_fun(dom,mu,v,n_ray,r,fun_level),...
            fun_span*[-1 1],fun_resol),num2cell(n,1),'un',0);
        init_sign=cell2mat(init_sign);
    end
end

end