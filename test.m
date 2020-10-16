mu=zeros(4,1);
v=eye(4);

fun=@(x1,x2,x3,x4) x1.^2+2*x2.^2+x3.^2+x4.^2-1;
mc_samples=round(10.^linspace(1,4,20));
plist=nan(size(mc_samples));
for i=1:length(mc_samples)
    i
    p=integrate_normal(mu,v,fun,'reg_type','fun','fun_span',3,'fun_resol',10,'mc_samples',mc_samples(i));
    plist(i)=p;
end

plot(mc_samples,plist,'-o')
set(gca,'xscale','log')

z=linspace(-10,10,100);
f1=(1-sign(z).*chi2cdf(z.^2,1))/2;
f2=normcdf(z,'upper');
plot(z,f1); hold on
plot(z,f2)

%%
load camouflage_edge_data
absent=[edge_powers_2(:,1),edge_lpr_2(:,1)];
present=[edge_powers_2(:,2),edge_lpr_2(:,2)];
n_samp=size(absent,1);
n_bootstrap=10;

absent_bootstrap=mvnrnd(mean(absent),cov(absent),n_samp*n_bootstrap);
% absent_bootstrap=reshape(absent_bootstrap,[n_samp,2,n_bootstrap]);

present_bootstrap=mvnrnd(mean(present),cov(present),n_samp*n_bootstrap);
% present_bootstrap=reshape(present_bootstrap,[n_samp,2,n_bootstrap]);

results=classify_normals(absent,present,'type','samp','samp_opt',false);
hold on; xlim([-1 4]); ylim([-2000 5000])
bd=results.norm_bd;

bd_shift=bd;
q0_family=linspace(200,400,40)';

samp_p12=nan(size(q0_family));
samp_p21=nan(size(q0_family));
norm_p12=nan(size(q0_family));
norm_p21=nan(size(q0_family));
norm_p12_bootstrap=nan(length(q0_family),n_bootstrap);
norm_p21_bootstrap=nan(length(q0_family),n_bootstrap);

for i=1:length(q0_family)
    bd_shift.q0=q0_family(i)
    plot_boundary(bd_shift,2,'plot_type','line','line_color',[0 0 1]);
    results_shift=classify_normals(absent,present,'type','samp','reg',bd_shift,'samp_opt',false,'bPlot',false);
    
    samp_p12(i)=results_shift.samp_errmat(1,2)/n_samp;
    samp_p21(i)=results_shift.samp_errmat(2,1)/n_samp;
    %     samp_p(i,:)=[q0_family(i) samp_p12 samp_p21];
    
    norm_p12(i)=results_shift.norm_errmat(1,2)/.5;
    norm_p21(i)=results_shift.norm_errmat(2,1)/.5;
    %     norm_p(i,:)=[q0_family(i) norm_p12 norm_p21];
    for j=1:n_bootstrap
        results_shift_bootstrap=classify_normals(absent_bootstrap((j-1)*n_samp+1:j*n_samp,:),present_bootstrap((j-1)*n_samp+1:j*n_samp,:),'type','samp','reg',bd_shift,'samp_opt',false,'bPlot',false);
        
        norm_p12_bootstrap(i,j)=results_shift_bootstrap.samp_errmat(1,2)/n_samp;
        norm_p21_bootstrap(i,j)=results_shift_bootstrap.samp_errmat(2,1)/n_samp;
    end
end

figure; hold on
plot(q0_family,samp_p12,'ob')
plot(q0_family,samp_p21,'or')
errorbar(q0_family,norm_p12,std(norm_p12_bootstrap,1,2),'-b')
errorbar(q0_family,norm_p21,std(norm_p21_bootstrap,1,2),'-r')

