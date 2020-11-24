function [p,pc,bd_pts]=int_norm_ray(mu,v,dom,varargin)
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
parser.KeepUnmatched=true;
addRequired(parser,'mu',@isnumeric);
addRequired(parser,'v',@isnumeric);
addRequired(parser,'dom');
addParameter(parser,'dom_type','quad');
addParameter(parser,'fun_span',5);
addParameter(parser,'fun_resol',100);
addParameter(parser,'AbsTol',1e-10);
addParameter(parser,'RelTol',1e-2);
addParameter(parser,'mc_samples',500);
parse(parser,mu,v,dom,varargin{:});
AbsTol=parser.Results.AbsTol;
RelTol=parser.Results.RelTol;
mc_samples=parser.Results.mc_samples;

dim=length(mu);

global bd_pts
bd_pts=[];

if dim==1
    p=int_norm_along_angles(mu,v,dom,varargin{:});
    pc=int_norm_along_angles(mu,v,dom,'side','complement',varargin{:});
elseif dim==2
    p=integral(@(theta) int_norm_along_angles(mu,v,dom,'theta',theta,varargin{:}),0,pi,'AbsTol',AbsTol,'RelTol',RelTol);
    pc=integral(@(theta) int_norm_along_angles(mu,v,dom,'theta',theta,'side','complement',varargin{:}),0,pi,'AbsTol',AbsTol,'RelTol',RelTol);
elseif dim==3
    p=integral2(@(theta,phi) int_norm_along_angles(mu,v,dom,'theta',theta,'phi',phi,varargin{:}),0,pi/2,0,2*pi,'AbsTol',AbsTol,'RelTol',RelTol);
    pc=integral2(@(theta,phi) int_norm_along_angles(mu,v,dom,'theta',theta,'phi',phi,'side','complement',varargin{:}),0,pi/2,0,2*pi,'AbsTol',AbsTol,'RelTol',RelTol);
elseif dim>3 % Monte-Carlo integration
    
    % uniform random rays (points on n-sphere)
    n_z=mvnrnd(zeros(dim,1),eye(dim),mc_samples)';
    n_z=n_z./vecnorm(n_z,2);
    
    p_rays=int_norm_along_rays(mu,v,dom,n_z,varargin{:})/2;
    p=mean(p_rays);
    pc_rays=int_norm_along_rays(mu,v,dom,n_z,'side','complement',varargin{:});
    pc=mean(pc_rays)/2;
end
end