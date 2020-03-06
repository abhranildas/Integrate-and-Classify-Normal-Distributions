mu_a=[0;0];
v=eye(2);

x=10.^linspace(0,3,100);
d=zeros(size(x));
for i=1:length(x)
    i
    mu_b=[x(i);0];
    results=classify_normals([mu_a,v],[mu_b,v],'bPlot',false);
    d(i)=results.d_norm;
end

plot(x,d)

%%

mu_a=[0;0]; mu_b=[100;0];
v=eye(2);

results=classify_normals([mu_a,v],[mu_b,1.01*v],'bPlot',false)