function [p_ray,bd_pts]=prob_theta(mu,v,reg_fn,theta,phi,side)

dim=length(mu);

if dim==1
    n_z=1;
elseif dim==2
    n_z=[cos(theta);sin(theta)];
elseif dim==3
    x=sin(theta).*cos(phi);
    y=sin(theta).*sin(phi);
    z=cos(theta);
    n_z=[x(:)';y(:)';z(:)'];
end

% unit rays in the original space:
n_x=sqrtm(v)*n_z;
n_x=n_x./vecnorm(n_x);

% initial signs and boundary distances along each direction
[init_sign,x]=reg_fn(n_x);

% relative boundary points in original space
bd_pts_rel=cellfun(@(x_ray,n_ray) x_ray.*n_ray, x,num2cell(n_x,1),'un',0);

if nargout==2
    % boundary points in original space
    bd_pts=horzcat(bd_pts_rel{:})+mu;
end

% standard boundary points
bd_pts_std=cellfun(@(bd_pt_rel) sqrtm(v)\bd_pt_rel, bd_pts_rel,'un',0);

% standard boundary distances, sorted
z=cellfun(@(n_ray,bd_pt_std) sort(n_ray'*bd_pt_std), num2cell(n_z,1), bd_pts_std,'un',0);

% total integral
if ~exist('side','var')
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
