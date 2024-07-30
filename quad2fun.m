function f=quad2fun(quad,fun_level,bPlot)
if nargin<3
    bPlot=0;
end

if ~bPlot
    % fast function of single vector argument
    f=@(x) dot(x,quad.q2*x) + quad.q1'*x + quad.q0 - fun_level;
else
    % slower function of multiple variables, for implicit plots
    dim=length(quad.q1);    
    if dim==1
        f=@(x) quad.q2*x.^2 + quad.q1*x + quad.q0 - fun_level;
    elseif dim==2
        f_each=@(x,y) [x y]*quad.q2*[x;y] + quad.q1'*[x;y] + quad.q0 - fun_level;
        f=@(x,y) arrayfun(f_each,x,y);
    elseif dim==3
        f_each=@(x,y,z) [x y z]*quad.q2*[x;y;z] + quad.q1'*[x;y;z] + quad.q0 - fun_level;
        f=@(x,y,z) arrayfun(f_each,x,y,z);
    end
end