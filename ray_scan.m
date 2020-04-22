function [init_sign,x,samp_correct]=ray_scan(reg,n,orig,varargin)

% parse inputs
parser = inputParser;

addRequired(parser,'reg',@(x) isstruct(x)|| isa(x,'function_handle'));
addRequired(parser,'n',@isnumeric);
addRequired(parser,'orig',@isnumeric);
addParameter(parser,'reg_type','quad');
addParameter(parser,'cheb_reg_span',3);
addParameter(parser,'func_crossings',100);

parse(parser,reg,n,orig,varargin{:});

reg_type=parser.Results.reg_type;
cheb_reg_span=parser.Results.cheb_reg_span;
func_crossings=parser.Results.func_crossings;

% pre-allocate variables in static workspace
% r=[];
% theta=[];
% az=[];
% el=[];
% 
% syms r theta az el

% root(s) along a ray through quad region
    function x_ray=roots_ray_quad(q2_pt,q1_pt)
        x_ray=sort(roots([q2_pt q1_pt a0]))';
        x_ray=x_ray(~imag(x_ray)); % only real roots
        
        % remove any roots that are tangents. Only crossing points
        slope=2*q2_pt*x_ray+q1_pt;
        x_ray=x_ray(slope~=0);
    end

% root(s) along a ray through cheb region
    function [init_sign_ray,r_ray]=roots_ray_cheb(n_ray)
        if nargin(reg)==1
            %r_ray=roots(reg(orig+n_ray*r))';
            %r_ray=roots(reg1)';
            [init_sign_ray,r_ray]=all_roots(reg_polar,cheb_reg_span*[-1 1],func_crossings);
            %syms r
            
%             if isempty(r_ray) % if no roots
%                 r_sign=0; 
%                 %init_sign_ray=sign(reg_polar(0)); % consider sign at 0
%             else
%                 r_sign=min(r_ray)-1e-1; 
%                 % consider sign of derivative at lowest root
%                 %init_sign_ray=-sign(reg_polar_diff(min(r_ray)));
%             end
%             init_sign_ray=sign(reg_polar(r_sign));
           % init_sign_ray=sign(reg(orig+n_ray*r_sign));
        elseif nargin(reg)==2
            %r_ray=roots(reg(orig(1)+n_ray(1)*r,orig(2)+n_ray(2)*r))';
            theta=cart2pol(n_ray(1),n_ray(2));
            %r_ray=roots(reg2(:,theta))';
            [init_sign_ray,r_ray]=all_roots(@(r) reg_polar(r,theta),cheb_reg_span*[-1 1],func_crossings);
%             if isempty(r_ray) % if no roots
%                 r_sign=0; % consider sign at 0
%                 %init_sign_ray=sign(reg_polar(0,theta));
%             else
%                 %init_sign_ray=-sign(reg_polar_diff(min(r_ray),theta));
%                 r_sign=min(r_ray)-1; % consider sign just south of lowest root
%             end
%             init_sign_ray=sign(reg_polar(r_sign,theta));
            %init_sign_ray=sign(reg(orig(1)+n_ray(1)*r_sign,orig(2)+n_ray(2)*r_sign));
        elseif nargin(reg)==3
            %r_ray=roots(reg(orig(1)+n_ray(1)*r,orig(2)+n_ray(2)*r,orig(3)+n_ray(3)*r))';
            [az,el]=cart2sph(n_ray(1),n_ray(2),n_ray(3));
            [init_sign_ray,r_ray]=all_roots(@(r) reg_polar(r,az,el),cheb_reg_span*[-1 1],func_crossings);
            %r_ray=roots(reg3(:,az,el))';
%             if isempty(r_ray) % if no roots
%                 r_sign=0; % consider sign at 0
%                 %init_sign_ray=sign(reg_polar(0,az,el));
%             else
%                 %init_sign_ray=-sign(reg_polar_diff(min(r_ray),az,el));
%                 r_sign=min(r_ray)-1; % consider sign just south of lowest root
%             end
            %init_sign_ray=sign(reg(orig(1)+n_ray(1)*r_sign,orig(2)+n_ray(2)*r_sign,orig(3)+n_ray(3)*r_sign));
%             init_sign_ray=sign(reg_polar(r_sign,az,el));
        end
    end

if strcmp(reg_type,'quad')
    
    if nargout==3
        [~,~,samp_correct]=samp_value(n',n',reg);
        init_sign=[];
        x=[];
    else
        n=n./vecnorm(n); % normalize direction vectors
        
        % boundary coefficients wrt origin
        a2=reg.a2;
        a1=2*reg.a2*orig+reg.a1;
        a0=orig'*reg.a2*orig+reg.a1'*orig+reg.a0;
        
        q2=dot(n,a2*n);
        q1=a1'*n;
        
        % sign of the quadratic at -inf:
        init_sign=sign(q2); % square term sets the sign
        init_sign(~init_sign)=-sign(q1(~init_sign)); % linear term sets the sign for leftovers
        init_sign(~init_sign)=sign(a0);% constant term sets the sign for the leftovers
        
        x=arrayfun(@roots_ray_quad,q2,q1,'un',0); % this allows function to calculate on multiple directions at once
    end
    
elseif strcmp(reg_type,'cheb')
    if nargout==3
        [~,~,samp_correct]=samp_value(n',n',reg,'reg_type','cheb');
        init_sign=[];
        x=[];
    else
        %syms r
        if nargin(reg)==1
            %r=chebfun('r',cheb_reg_span*[-1 1],'splitting','on');
            reg_polar=@(r) reg(orig+r);            
            %reg_polar_diff=matlabFunction(diff(reg_polar(r)));
            %reg1=chebfun(@(r) reg(orig+r),cheb_reg_span*[-1 1]);
        elseif nargin(reg)==2
            
            %reg2=chebfun2(reg,repelem(orig',2)+cheb_reg_span*[-1 1 -1 1]);
            %reg2_polar=chebfun2(@(r,theta) reg2(orig(1)+r.*cos(theta),orig(2)+r.*sin(theta)));
            %reg2=chebfun2(@(r,theta) reg(orig(1)+r.*cos(theta),orig(2)+r.*sin(theta)),[cheb_reg_span*[-1 1] 0 pi]);
            reg_polar=@(r,theta) reg(orig(1)+r.*cos(theta),orig(2)+r.*sin(theta));
            %reg_polar_diff=matlabFunction(diff(reg_polar(r,theta)));
            %             reg2c=chebfun2(@(x,y) reg(orig(1)+x,orig(2)+y),cheb_reg_span*[-1 1 -1 1]);
        elseif nargin(reg)==3
            %reg3=chebfun3(@(r,az,el) reg(orig(1)+r.*cos(el).*cos(az),orig(2)+r.*cos(el).*sin(az),orig(3)+r.*sin(el)),[cheb_reg_span*[-1 1] 0 2*pi 0 pi/2]);
            reg_polar=@(r,az,el) reg(orig(1)+r.*cos(el).*cos(az),orig(2)+r.*cos(el).*sin(az),orig(3)+r.*sin(el));
        %reg_polar_diff=matlabFunction(diff(reg_polar(r,az,el),r));
        end
        [init_sign,x]=cellfun(@roots_ray_cheb,num2cell(n,1),'un',0); % this allows function to calculate on multiple directions at once
        init_sign=cell2mat(init_sign);
    end
end

end