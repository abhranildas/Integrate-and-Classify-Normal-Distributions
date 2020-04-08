function points=fibonacci_sphere(n_pts)
% Fibonacci sphere algorithm for evenly distributed
% points on a sphere, translated from:
% https://stackoverflow.com/a/26127012/711017
% Credits:
%   Abhranil Das <abhranil.das@utexas.edu>
%	R Calen Walshe
%	Wilson S Geisler
%	Center for Perceptual Systems, University of Texas at Austin
% If you use this code, please cite:
%   A new method to compute classification error
%   https://jov.arvojournals.org/article.aspx?articleid=2750251

points=nan(3,n_pts);
offset=2./n_pts;
increment=pi*(3-sqrt(5));
for i=0:n_pts-1
    y=i*offset-1 + offset/2;
    r=sqrt(1-y^2);
    phi=mod(i+1,n_pts)*increment;
    x=r*cos(phi);
    z=r*sin(phi);
    points(:,i+1)=[x;y;z];
end