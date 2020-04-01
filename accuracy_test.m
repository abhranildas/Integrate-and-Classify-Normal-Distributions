%% Classification converging with n_rays
mu_1=[0;0;0];
v_1=eye(3);

mu_2=[2;1;1];
v_2=2*eye(3);

results=classify_normals([mu_1,v_1],[mu_2,v_2],'bPlot',false);
gx2_estimate=results.norm_err;

% force ray method by supplying boundary functions
reg_fn_1=@(n) opt_reg(n,[mu_1,v_1],[mu_2,v_2]);  % optimum boundary function for normal 1
reg_fn_2=@(n) opt_reg(n,[mu_2,v_2],[mu_1,v_1]);  % optimum boundary function for normal 2

n_rays=round(10.^linspace(1,5,15));
p_ray=nan(size(n_rays));
for i=1:length(n_rays)
results_ray=classify_normals([mu_1,v_1],[mu_2,v_2],'custom_reg',{reg_fn_1,reg_fn_2},'n_rays',n_rays(i),'bPlot',false);
p_ray(i)=results_ray.norm_err
end
rel_errs=(p_ray-gx2_estimate)/gx2_estimate;
plot(n_rays,abs(rel_errs),'-o')
set(gca,'xscale','log')
set(gca,'yscale','log')

%% Integration converging with n_rays
% TODO second x-axis, just below the first, of process time (since prop. to n_rays)
mu=[0;0;0];
v=eye(3);

reg_coeffs.a2=-0.5*eye(3);
reg_coeffs.a1=[-2;-1;-1];
reg_coeffs.a0=5;

p_gx2=integrate_normal(mu,v,reg_coeffs,'bPlot',false);

% force ray method by supplying boundary functions
reg_fn=@(n) quad_reg(n,[mu,v],reg_coeffs);  % optimum boundary function for normal 1

n_rays=round(10.^linspace(1,7,30));
p_ray=nan(size(n_rays));
time=nan(size(n_rays));
for i=1:length(n_rays)
    p_ray(i)=integrate_normal(mu,v,reg_fn,'n_rays',n_rays(i),'bPlot',false);
    time(i)=timeit(@() integrate_normal(mu,v,reg_fn,'n_rays',n_rays(i),'bPlot',false));
    i
end

rel_errs=(p_ray-p_gx2)/p_gx2;

yyaxis left
plot(n_rays,abs(rel_errs),'-o')
set(gca,'xscale','log')
set(gca,'yscale','log')
ylabel 'relative error'

yyaxis right
plot(n_rays,time,'-o')
ylabel 'process time (s)'

xlabel 'n_{rays}'

%% performance vs distance

mu_1=[0;0;0];
v=eye(3);

d_true=10.^linspace(0,3,500);
d_gx2=zeros(size(d_true));
d_ray=zeros(size(d_true));
for i=1:length(d_true)
    i
    mu_2=[d_true(i);0;0];
    results_gx2=classify_normals([mu_1,v],[mu_2,v],'n_rays',1,'bPlot',false);
    d_gx2(i)=results_gx2.norm_d;
    
    % force ray method by supplying boundary functions
    reg_fn_1=@(n) opt_reg(n,[mu_1,v],[mu_2,v]);  % optimum boundary function for normal 1
    reg_fn_2=@(n) opt_reg(n,[mu_2,v],[mu_1,v]);  % optimum boundary function for normal 2
    results_ray=classify_normals([mu_1,v],[mu_2,v],'custom_reg',{reg_fn_1,reg_fn_2},'n_rays',1e4,'bPlot',false);
    d_ray(i)=results_ray.norm_d;
end

rel_err_gx2=(d_gx2-d_true)./d_true;
rel_err_ray=(d_ray-d_true)./d_true;

figure
plot(d_true,abs(rel_err_gx2),'-o')
hold on
plot(d_true,abs(rel_err_ray),'-o')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel 'true d'''
ylabel 'relative error'

%%

mu_1=[0;0]; mu_2=[100;0];
v=eye(2);

results=classify_normals([mu_1,v],[mu_2,1.01*v],'bPlot',false)