function [p,pc,bd_pts]=int_norm_grid(mu,v,reg_fn,n_rays)
% Integrate a normal distribution (upto 3D) over a specified region, using
% the ray method.
%
% How to use this command:
% See github readme at https://github.com/abhranildas/classify
%
% Credits:
%   Abhranil Das <abhranil.das@utexas.edu>
%	Wilson S Geisler
%	Center for Perceptual Systems, University of Texas at Austin
% If you use this code, please cite:
%   A new method to compute classification error
%   https://jov.arvojournals.org/article.aspx?articleid=2750251

dim=length(mu);

if ~exist('n_rays','var')
    n_rays=1e4;
end

% Create grid of unit vectors
if dim==1
    n_z=1;
elseif dim==2
    dth=pi/n_rays;
    th=0:dth:pi;
    n_z=[cos(th);sin(th)];
elseif dim==3
    n_z=fibonacci_sphere(n_rays);
    %     dth=5e-2; dph=5e-2;
    %     th=0:dth:pi/2;
    %     ph=0:dph:2*pi;
    %     n_list=[];
    %     for i=1:length(th)
    %         for j=1:length(ph)
    %             n_list=[n_list,[sin(th(i))*cos(ph(j)); sin(th(i))*sin(ph(j)); cos(th(i))]];
    %         end
    %     end
end

% unit rays in the original space:
n_x=sqrtm(v)*n_z;
n_x=n_x./vecnorm(n_x);

% initial signs and boundary distances along each direction
[init_sign,x]=reg_fn(n_x);

% relative boundary points in original space
bd_pts_rel=cellfun(@(x_ray,n_ray) x_ray.*n_ray, x,num2cell(n_x,1),'un',0);

% boundary points in original space
bd_pts=horzcat(bd_pts_rel{:})+mu;

% standard boundary points
bd_pts_std=cellfun(@(bd_pt_rel) sqrtm(v)\bd_pt_rel, bd_pts_rel,'un',0);

% standard boundary distances, sorted
z=cellfun(@(n_ray,bd_pt_std) sort(n_ray'*bd_pt_std), num2cell(n_z,1), bd_pts_std,'un',0);

    % function to compute probability slice along each ray
    function [p_ray,p_ray_c]=prob_ray(init_sign_ray,z_ray)
        [f_big,f_small]=norm_rad_surv_split(z_ray,dim);
        p_ray_big=(init_sign_ray+1)/2+init_sign_ray*sum((-1).^(1:length(z_ray)).*f_big);
        p_ray_small=init_sign_ray*sum((-1).^(1:length(z_ray)).*f_small);
        p_ray=p_ray_big+p_ray_small;
        p_ray_c=1-p_ray_big-p_ray_small;
        %norm_rad_surv=(1-sign(z_ray).*chi2cdf(z_ray.^2,dim))/2; % normal radial survival function
        %p_ray=1/2+init_sign_ray*(1/2+ sum((-1).^(1:length(z_ray)).*norm_rad_surv));
    end

% total integral
[p_rays,pc_rays]=cellfun(@prob_ray,num2cell(init_sign),z);
p=mean(p_rays);
pc=mean(pc_rays);
end