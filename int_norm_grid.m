function [mass,bd_pts]=int_norm_grid(mu,v,bd_fn,n_points)
% Integrate a normal distribution (upto 3D) over a specified region, using
% the grid method.
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

% Create grid of unit vectors
if dim==1
    n_list=[1];
elseif dim==2
    dth=pi/n_points;
    th=0:dth:pi;
    n_list=[cos(th);sin(th)];    
elseif dim==3
    n_list=fibonacci_sphere(n_points);
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

% Integrate
bd_pts_std=[];
bd_sign=[];

[~,mass]=bd_fn(zeros(dim,1)); % starting mass =1 if mu is in region.

m=0;
[r,r_sign]=bd_fn(n_list);
bd_pts_std=cellfun(@(x,y) x.*y,r,num2cell(n_list,1),'un',0);
bd_pts_std=horzcat(bd_pts_std{:});
r=horzcat(r{:});
r_sign=horzcat(r_sign{:});
if dim==1
    dm=sum(r_sign.*normcdf(abs(r),'upper'));
elseif dim==2
    dm=sum(r_sign.*exp(-r.^2/2));
elseif dim==3
    %sinth=sqrt(1-n(3)^2);
    dm=sum(r_sign.*(sqrt(pi)/2 * (1-erf(abs(r)/sqrt(2))) + abs(r).*exp(-abs(r).^2/2)/sqrt(2)));%*sinth;
end
m=m+dm;
% for i=1:length(n_list)
%     n=n_list(:,i);
%     [r,r_sign]=bd_fn(n);
%         bd_pt_std=r.*n; % standardized boundary point
%         bd_pts_std=[bd_pts_std,bd_pt_std];
%         bd_sign=[bd_sign,r_sign];
%         if dim==1            
%             dm=sum(r_sign.*normcdf(abs(r),'upper'));
%         elseif dim==2
%             dm=sum(r_sign.*exp(-r.^2/2));
%         elseif dim==3
%             %sinth=sqrt(1-n(3)^2);
%             dm=sum(r_sign.*(sqrt(pi)/2 * (1-erf(abs(r)/sqrt(2))) + abs(r).*exp(-abs(r).^2/2)/sqrt(2)));%*sinth;
%         end
%         m=m+dm;
% end

% factors
if dim==2
    m=m*dth/(2*pi);
elseif dim==3
    %m=m*dth*dph/(2*pi^(3/2));
    m=m/(n_points*sqrt(pi));
end
mass=mass+m;

% un-standardize boundary points:
C=chol(v,'lower'); % Cholesky decomposition of vcov
if numel(bd_pts_std)
    bd_pts=C*bd_pts_std+mu;
end