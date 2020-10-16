function [init_sign,r]=ray_scan_fun(fun,xrange,resolution)
% Inputs :
% f : function of one variable
% [xmin - xmax] : range where f is continuous containing zeros
% N : control of the minimum distance (xmax-xmin)/N between two zeros

% g=sym(f); %converted to symbolic
% h=finverse(g); %inverse of the function
% F=matlabFunction(h)

if (nargin<3)
    resolution=1e2;
end

x=linspace(xrange(1),xrange(2),resolution);
y=fun(x);
init_sign=sign(y(1));
sign_diff=diff(sign(y));

sign_change_locs=find(sign_diff);
if isempty(sign_change_locs)    
    r=double.empty(1,0);
else
    n_sign_change=numel(sign_change_locs);
    
    r=nan(1,n_sign_change);
    for i=1:n_sign_change
        r(i)=fzero(fun,x(sign_change_locs(i)+[0 1]));
    end
end
r=unique(r);