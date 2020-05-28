%% Generalized chi-squared demo
% look into each function code for more documentation

lambda=[1 10 2];
m=[1 2 3];
delta=[2 3 7];
c=10;

% calculate PDF and CDF at a point
x=25;
format long
f=gx2pdf(x,lambda,m,delta,c)
p_imhof=gx2cdf_imhof(x,lambda,m,delta,c) % Imhof's method (recommended)
p_ruben=gx2cdf_ruben(x,lambda,m,delta,c) % Ruben's method

% plot PDF and CDF
lambda=[1 -10 2]; % mixed coefficients (Ruben's method cannot handle this)
[mu,v]=gx2stat(lambda,m,delta,c) % mean and variance
x=linspace(mu-3*sqrt(v),mu+3*sqrt(v),1e3);
f=arrayfun(@(x) gx2pdf(x,lambda,m,delta,c),x);
p=arrayfun(@(x) gx2cdf_imhof(x,lambda,m,delta,c),x);

yyaxis left
area(x,f,'facealpha',.5)
xline(mu,'-',{'mean'});
ylim([0 .015])
ylabel 'PDF'
yyaxis right
area(x,p,'facealpha',.5)
ylabel 'CDF'
xlabel x
xlim([mu-3*sqrt(v),mu+3*sqrt(v)])