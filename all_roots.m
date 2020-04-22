function [init_sign,r]=all_roots(f,xrange,f_crossings)
% Inputs :
% f : function of one variable
% [xmin - xmax] : range where f is continuous containing zeros
% N : control of the minimum distance (xmax-xmin)/N between two zeros
if (nargin<3)
    f_crossings=1e2;
end

x=linspace(xrange(1),xrange(2),f_crossings);
y=f(x);
sign_diff=diff(sign(y));
%first_root_loc=
[~,~,init_diff]=find(sign_diff,1);
if isempty(init_diff)
    init_sign=sign(y(1));
    r=double.empty(1,0);
else
    init_sign=-sign(init_diff);
    sign_change_locs=find(sign_diff);
    n_sign_change=numel(sign_change_locs);
    
    r=nan(1,n_sign_change);
    for i=1:n_sign_change
        r(i)=fzero(f,x(sign_change_locs(i)+[0 1]));
    end
end

% for i=1:n
%     x1=x2;
%     y1=y2;
%     x2=xmin+i*dx;
%     y2=f(x2);
%     if (y1*y2<=0)                              % Rolle's theorem : one zeros (or more) present
%         r(i)=fsolve(f,(x2*y1-x1*y2)/(y1-y2),optimoptions('fsolve','Display','none')); % Linear approximation to guess the initial value in the [x1,x2] range.
%     end
% end
%r=unique(r(~isnan(r)));
%r=unique(r);