function [init_sign,x]=ray_scan(reg,reg_type,n,orig)

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
            r_ray=roots(reg(orig+n_ray*r))';
            if isempty(r_ray) % if no roots
                r_sign=0; % consider sign at 0
            else
                r_sign=min(r_ray)-1; % consider sign just south of lowest root
            end
            init_sign_ray=sign(reg(orig+n_ray*r_sign));
        elseif nargin(reg)==2
            r_ray=roots(reg(orig(1)+n_ray(1)*r,orig(2)+n_ray(2)*r))';
            if isempty(r_ray) % if no roots
                r_sign=0; % consider sign at 0
            else
                r_sign=min(r_ray)-1; % consider sign just south of lowest root
            end
            init_sign_ray=sign(reg(orig(1)+n_ray(1)*r_sign,orig(2)+n_ray(2)*r_sign));
        elseif nargin(reg)==3
            r_ray=roots(reg(orig(1)+n_ray(1)*r,orig(2)+n_ray(2)*r,orig(3)+n_ray(3)*r))';
            if isempty(r_ray) % if no roots
                r_sign=0; % consider sign at 0
            else
                r_sign=min(r_ray)-1; % consider sign just south of lowest root
            end
            init_sign_ray=sign(reg(orig(1)+n_ray(1)*r_sign,orig(2)+n_ray(2)*r_sign,orig(3)+n_ray(3)*r_sign));
        end
    end

if strcmp(reg_type,'quad')
    
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
    
elseif strcmp(reg_type,'cheb')
    r=chebfun('r',[-inf inf],'splitting','on');
    [init_sign,x]=cellfun(@roots_ray_cheb,num2cell(n,1),'un',0); % this allows function to calculate on multiple directions at once
    init_sign=cell2mat(init_sign);
end

end