function [init_sign,x]=quad_reg(n,coeffs,orig)
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

% boundary coefficients wrt origin
a2=coeffs.a2;
a1=2*coeffs.a2*orig+coeffs.a1;
a0=orig'*coeffs.a2*orig+coeffs.a1'*orig+coeffs.a0;

n=n./vecnorm(n); % normalize direction vectors
q2=dot(n,a2*n);
q1=a1'*n;

% sign of the quadratic at -inf:
init_sign=sign(q2); % square term sets the sign
init_sign(~init_sign)=-sign(q1(~init_sign)); % linear term sets the sign for leftovers
init_sign(~init_sign)=sign(a0);% constant term sets the sign for the leftovers

    % function to return root(s) along a ray
    function x_ray=roots_ray(q2_pt,q1_pt)
        x_ray=sort(roots([q2_pt q1_pt a0]))';
        x_ray=x_ray(~imag(x_ray)); % only real roots
        
        % remove any roots that are tangents. Only crossing points
        slope=2*q2_pt*x_ray+q1_pt;
        x_ray=x_ray(slope~=0);
    end

x=arrayfun(@roots_ray,q2,q1,'un',0); % this allows function to calculate on multiple directions at once
end