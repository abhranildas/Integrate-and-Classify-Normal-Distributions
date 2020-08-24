lambda=1;
m=2;
delta=0;
sigma=0;
c=0;

x=linspace(-1,10,100);
f1=arrayfun(@(x)gx2cdf_davies(x,lambda,m,delta,sigma,c),x);

fplot(@(x) gx2cdf(x,lambda,m,delta,sigma,c));

fplot(@(x) gx2pdf(x,lambda,m,delta,sigma,c));

f1=arrayfun(@(x)ncx2pdf(x,lambda,m,delta),x);
f2=arrayfun(@(x)gx2pdf(x,lambda,m,delta,sigma,c),x);

hold on
plot(x,f1)
plot(x,f2)
plot(x,f3)

[p,errbnd]=gx2cdf_ruben(x,lambda,m,delta,c,1000);
