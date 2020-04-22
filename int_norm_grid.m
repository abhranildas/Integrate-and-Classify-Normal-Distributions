function [p,pc,bd_pts]=int_norm_grid(mu,v,reg,varargin)
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

parser = inputParser;
addRequired(parser,'mu',@isnumeric);
addRequired(parser,'v',@isnumeric);
addRequired(parser,'reg');
addParameter(parser,'reg_type','quad');
addParameter(parser,'cheb_reg_span',3);
addParameter(parser,'func_crossings',100);
addParameter(parser,'AbsTol',1e-10);
addParameter(parser,'RelTol',1e-2);
parse(parser,mu,v,reg,varargin{:});
reg=parser.Results.reg;
reg_type=parser.Results.reg_type;
cheb_reg_span=parser.Results.cheb_reg_span;
func_crossings=parser.Results.func_crossings;
AbsTol=parser.Results.AbsTol;
RelTol=parser.Results.RelTol;

dim=length(mu);

% if ~exist('n_rays','var')
%     n_rays=1e4;
% end

% Grid of unit vectors in standardized space
% if dim==1
%     n_z=1;
% elseif dim==2
%     dth=pi/n_rays;
%     th=0:dth:pi;
%     n_z=[cos(th);sin(th)];
% elseif dim==3
%     n_z=fibonacci_sphere(n_rays);
%     %     dth=5e-2; dph=5e-2;
%     %     th=0:dth:pi/2;
%     %     ph=0:dph:2*pi;
%     %     n_list=[];
%     %     for i=1:length(th)
%     %         for j=1:length(ph)
%     %             n_list=[n_list,[sin(th(i))*cos(ph(j)); sin(th(i))*sin(ph(j)); cos(th(i))]];
%     %         end
%     %     end
% end

% unit rays in the original space:
% n_x=sqrtm(v)*n_z;
% n_x=n_x./vecnorm(n_x);

% initial signs and boundary distances along each direction
% [init_sign,x]=reg_fn(n_x);

% relative boundary points in original space
% bd_pts_rel=cellfun(@(x_ray,n_ray) x_ray.*n_ray, x,num2cell(n_x,1),'un',0);

% boundary points in original space
% bd_pts=horzcat(bd_pts_rel{:})+mu;

% standard boundary points
% bd_pts_std=cellfun(@(bd_pt_rel) sqrtm(v)\bd_pt_rel, bd_pts_rel,'un',0);

% standard boundary distances, sorted
% z=cellfun(@(n_ray,bd_pt_std) sort(n_ray'*bd_pt_std), num2cell(n_z,1), bd_pts_std,'un',0);

% total integral
% [p_rays,pc_rays]=cellfun(@(init_sign_ray,z_ray) prob_ray(init_sign_ray,z_ray,dim),num2cell(init_sign),z);
% p=mean(p_rays);
% pc=mean(pc_rays);

global bd_pts
bd_pts=[];

if dim==1
    p=prob_bd_angle(mu,v,reg,'reg_type',reg_type,'cheb_reg_span',cheb_reg_span,'func_crossings',func_crossings);
    pc=prob_bd_angle(mu,v,reg,'reg_type',reg_type,'cheb_reg_span',cheb_reg_span,'func_crossings',func_crossings,'side','complement');
elseif dim==2
    p=integral(@(theta) prob_bd_angle(mu,v,reg,'reg_type',reg_type,'cheb_reg_span',cheb_reg_span,'func_crossings',func_crossings,'theta',theta),0,pi,'AbsTol',AbsTol,'RelTol',RelTol);
    pc=integral(@(theta) prob_bd_angle(mu,v,reg,'reg_type',reg_type,'cheb_reg_span',cheb_reg_span,'func_crossings',func_crossings,'theta',theta,'side','complement'),0,pi);    
elseif dim==3
    p=integral2(@(theta,phi) prob_bd_angle(mu,v,reg,'reg_type',reg_type,'cheb_reg_span',cheb_reg_span,'func_crossings',func_crossings,'theta',theta,'phi',phi),0,pi/2,0,2*pi,'AbsTol',AbsTol,'RelTol',RelTol);
    pc=integral2(@(theta,phi) prob_bd_angle(mu,v,reg,'reg_type',reg_type,'cheb_reg_span',cheb_reg_span,'func_crossings',func_crossings,'theta',theta,'phi',phi,'side','complement'),0,pi/2,0,2*pi,'AbsTol',AbsTol,'RelTol',RelTol);
    %[~,bd_pts]=prob_bd_angle(mu,v,reg_fn,linspace(0,pi,1e4),nan);
end
% p=ps(1); pc=ps(2);
end