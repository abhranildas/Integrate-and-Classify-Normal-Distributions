function [init_sign,x,samp_correct]=ray_scan(reg,n,varargin)

% parse inputs
dim=size(n,1);
parser = inputParser;
addRequired(parser,'reg',@(x) isstruct(x)|| isa(x,'function_handle'));
addRequired(parser,'n',@isnumeric);
addParameter(parser,'mu',zeros(dim,1),@isnumeric);
addParameter(parser,'v',eye(dim),@isnumeric);
addParameter(parser,'reg_type','quad');
addParameter(parser,'fun_span',3);
addParameter(parser,'fun_resol',100);

parse(parser,reg,n,varargin{:});
mu=parser.Results.mu;
v=parser.Results.v;
reg_type=parser.Results.reg_type;
fun_span=parser.Results.fun_span;
fun_resol=parser.Results.fun_resol;

    % root(s) along a ray through quad region
    function x_ray=roots_ray_quad(q2_pt,q1_pt)
        x_ray=sort(roots([q2_pt q1_pt q0]))';
        x_ray=x_ray(~imag(x_ray)); % only real roots
        
        % remove any roots that are tangents. Only crossing points
        slope=2*q2_pt*x_ray+q1_pt;
        x_ray=x_ray(slope~=0);
    end

    % root(s) along a ray through fun region
    function [init_sign_ray,r_ray]=fun_roots_ray(n_ray,reg)
%         if dim==1
            [init_sign_ray,r_ray]=all_roots(@(r) standard_ray_fun(reg,mu,v,n_ray,r),fun_span*[-1 1],fun_resol);
%         elseif dim==2
            %theta=cart2pol(n_ray(1),n_ray(2));
%             [init_sign_ray,r_ray]=all_roots(@(r) standard_ray_func(reg,mu,v,n_ray,r),func_span*[-1 1],func_crossings);
%         elseif dim==3
            %[az,el]=cart2sph(n_ray(1),n_ray(2),n_ray(3));
%             [init_sign_ray,r_ray]=all_roots(@(r) standard_ray_func(reg,mu,v,n_ray,r),func_span*[-1 1],func_crossings);
%         end
    end


if strcmp(reg_type,'quad')
    
    if nargout==3
        [~,~,samp_correct]=samp_value(n',n',reg);
        init_sign=[];
        x=[];
    else
        n=n./vecnorm(n); % normalize direction vectors
        
        % standardized boundary coefficients
        quad_s=standard_quad(reg,mu,v);
        %q2=reg.q2;
        %q1=2*reg.q2*mu+reg.q1;
        %q0=mu'*reg.q2*mu+reg.q1'*mu+reg.q0;
        
        q2=dot(n,quad_s.q2*n);
        q1=quad_s.q1'*n;
        q0=quad_s.q0;
        
        % sign of the quadratic at -inf:
        init_sign=sign(q2); % square term sets the sign
        init_sign(~init_sign)=-sign(q1(~init_sign)); % linear term sets the sign for leftovers
        init_sign(~init_sign)=sign(q0);% constant term sets the sign for the leftovers
        
        x=arrayfun(@roots_ray_quad,q2,q1,'un',0); % this allows function to calculate on multiple directions at once
    end
    
elseif strcmp(reg_type,'fun')
    if nargout==3
        [~,~,samp_correct]=samp_value(n',n',reg,'reg_type','fun');
        init_sign=[];
        x=[];
    else
        [init_sign,x]=cellfun(@(n_ray) fun_roots_ray(n_ray,reg),num2cell(n,1),'un',0);
        init_sign=cell2mat(init_sign);
    end
end

end