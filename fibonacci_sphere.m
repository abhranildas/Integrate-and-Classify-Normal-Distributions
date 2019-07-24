function points=fibonacci_sphere(n_samples)
% Fibonacci sphere algorithm for evenly distributed
% points on a sphere, translated from:
% https://stackoverflow.com/a/26127012/711017
% to MATLAB by: Abhranil Das <abhranil.das@utexas.edu>
% Please cite if you use this code.

points=nan(3,n_samples);
offset=2./n_samples;
increment=pi*(3-sqrt(5));
for i=0:n_samples-1
    y=i*offset-1 + offset/2;
    r=sqrt(1-y^2);
    phi=mod(i+1,n_samples)*increment;
    x=r*cos(phi);
    z=r*sin(phi);
    points(:,i+1)=[x;y;z];
end