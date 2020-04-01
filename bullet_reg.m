function [init_sign,x]=bullet_reg(n,e,orig)
% Return distances and entry/exit signs of a quadratic region boundary (defined by
% coeffs), standardized wrt norm, in the direction(s) of ray(s) n.
% If n doesn't hit the boundary, sign indicates if this direction is inside
% or outside the boundary.
% For n=0, sign indicates whether the origin (mu) is inside, outside, or on the boundary.
%
% How to use this command:
% See github readme at https://github.com/abhranildas/classify
%
% Author:
%   Abhranil Das <abhranil.das@utexas.edu>
%	Center for Perceptual Systems, University of Texas at Austin
% If you use this code, please cite:
%   A new method to compute classification error
%   https://jov.arvojournals.org/article.aspx?articleid=2750251

n=n./vecnorm(n); % normalize direction vectors
ox=orig(1); oy=orig(2);

    % function to return root(s) along a ray
    function x_ray=roots_ray(n_ray)
        nx=n_ray(1); ny=n_ray(2);
        x_ray=roots([ny^2*nx, ny^2*ox+2*nx*ny*oy, nx*oy^2+2*ny*ox*oy, ox*oy^2-2*e])';
        x_ray=x_ray(~imag(x_ray)); % only real roots
    end

x=cellfun(@roots_ray,num2cell(n,1),'un',0); % this allows function to calculate on multiple directions at once
init_sign=-sign(n(1,:));
end