function r=all_roots0(f,xmin,xmax,n)
if ~exist('n','var')
    n=100;
end
r=nan(n,1);             % Initializes x.
starting_points=linspace(xmin,xmax,n);
for i=1:n
    % Look for the zeros in the function's current window.
    r(i)=fzero(f, starting_points(i));
end
r=r(diff(r)>1e-12);