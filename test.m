%%
load camouflage_edge_data
absent=[edge_powers_2(:,1),edge_lpr_2(:,1)];
present=[edge_powers_2(:,2),edge_lpr_2(:,2)];
n_samp=size(absent,1);
n_bootstrap=10;

% absent_bootstrap=mvnrnd(mean(absent),cov(absent),n_samp*n_bootstrap);
% absent_bootstrap=reshape(absent_bootstrap,[n_samp,2,n_bootstrap]);

% present_bootstrap=mvnrnd(mean(present),cov(present),n_samp*n_bootstrap);
% present_bootstrap=reshape(present_bootstrap,[n_samp,2,n_bootstrap]);

results=classify_normals(absent,present,'type','samp','samp_opt',false);
hold on; xlim([-1 4]); ylim([-2000 5000])
bd=results.norm_bd;

bd_shift=bd;
q0s=linspace(200,400,40)';

samp_p11=nan(size(q0s));
samp_p22=nan(size(q0s));
norm_p11=nan(size(q0s));
norm_p22=nan(size(q0s));
% norm_p12_bootstrap=nan(length(q0_family),n_bootstrap);
% norm_p21_bootstrap=nan(length(q0_family),n_bootstrap);

for i=1:length(q0s)
    bd_shift.q0=q0s(i)
    plot_boundary(bd_shift,2,'plot_type','line','line_color',[0 0 1]);
    results_shift=classify_normals(absent,present,'type','samp','reg',bd_shift,'samp_opt',false,'bPlot',false);
    
    samp_p11(i)=results_shift.samp_errmat(1,2)/n_samp;
    samp_p22(i)=results_shift.samp_errmat(2,1)/n_samp;
    %     samp_p(i,:)=[q0_family(i) samp_p12 samp_p21];
    
    norm_p11(i)=results_shift.norm_errmat(1,2)/.5;
    norm_p22(i)=results_shift.norm_errmat(2,1)/.5;
    
    %     norm_p(i,:)=[q0_family(i) norm_p12 norm_p21];
%     for j=1:n_bootstrap
%         results_shift_bootstrap=classify_normals(absent_bootstrap((j-1)*n_samp+1:j*n_samp,:),present_bootstrap((j-1)*n_samp+1:j*n_samp,:),'type','samp','reg',bd_shift,'samp_opt',false,'bPlot',false);
%         
%         norm_p12_bootstrap(i,j)=results_shift_bootstrap.samp_errmat(1,2)/n_samp;
%         norm_p21_bootstrap(i,j)=results_shift_bootstrap.samp_errmat(2,1)/n_samp;
%     end
end
norm_p11_sd=sqrt(n_samp*norm_p11.*(1-norm_p11))/n_samp;
norm_p22_sd=sqrt(n_samp*norm_p22.*(1-norm_p22))/n_samp;


figure; hold on
plot(q0s,samp_p11,'ob')
plot(q0s,samp_p22,'or')
errorbar(q0s,norm_p11,norm_p11_sd,'-b')
errorbar(q0s,norm_p22,norm_p22_sd,'-r')

%%
load camouflage_edge_data
absent=[edge_powers_2(:,1),edge_lpr_2(:,1),edge_powers_4(:,1),edge_lpr_4(:,1),edge_powers_8(:,1),edge_lpr_8(:,1)];
present=[edge_powers_2(:,2),edge_lpr_2(:,2),edge_powers_4(:,2),edge_lpr_4(:,2),edge_powers_8(:,2),edge_lpr_8(:,2)];
n_samp=1e2;
absent=absent(1:n_samp,:);
present=present(1:n_samp,:);

results=classify_normals(absent,present,'type','samp','samp_opt',false,'bPlot',false);
bd=results.norm_bd;


% shift criterion
bd_shift=bd;
q0s=linspace(300,475,50)';

samp_p11=nan(size(q0s));
samp_p22=nan(size(q0s));
norm_p11=nan(size(q0s));
norm_p22=nan(size(q0s));

for i=1:length(q0s)
    bd_shift.q0=q0s(i);
    results_shift=classify_normals(absent,present,'type','samp','reg',bd_shift,'samp_opt',false,'bPlot',false);
    
    samp_p11(i)=results_shift.samp_errmat(1,1)/n_samp;
    samp_p22(i)=results_shift.samp_errmat(2,2)/n_samp;
    
    norm_p11(i)=results_shift.norm_errmat(1,1)/.5;
    norm_p22(i)=results_shift.norm_errmat(2,2)/.5;    
end
norm_p11_sd=sqrt(norm_p11.*(1-norm_p11)/n_samp);
norm_p22_sd=sqrt(norm_p22.*(1-norm_p22)/n_samp);


figure; hold on
xline(bd.q0,'k','optimal criterion','LabelVerticalAlignment','bottom') % optimal criterion

% p11
x=q0s;
y=norm_p11;
dy=norm_p11_sd;
fill([x;flipud(x)],[y-dy;flipud(y+dy)],[0 0 1],'facealpha',.3,'linestyle','none');
plot(q0s,samp_p11,'.b','markersize',6)

% p22
y=norm_p22;
dy=norm_p22_sd;
fill([x;flipud(x)],[y-dy;flipud(y+dy)],[1 0 0],'facealpha',.3,'linestyle','none');
plot(q0s,samp_p22,'.r','markersize',6)

xlim([min(q0s) max(q0s)]); ylim([0 1])
xlabel('assumed criterion (q_0)')
set(gca,'fontsize',13); box off

% shift d'
dprime=results.norm_dprime;

mu_a=mean(absent)'; mu_b=mean(present)';
v_a=cov(absent); v_b=cov(present);
% center=(mu_a+mu_b)/2;
% delta=(mu_b-mu_a)/2;

v_scales=10.^linspace(-.5,2,50)';

dprimes=nan(size(v_scales));
samp_p11=nan(size(v_scales));
samp_p22=nan(size(v_scales));
norm_p11=nan(size(v_scales));
norm_p22=nan(size(v_scales));

for i=1:length(v_scales)
%     mu_a_shift=center-cov_scales(i)*delta;
%     mu_b_shift=center+cov_scales(i)*delta;
    results_pretend=classify_normals([mu_a,v_scales(i)*v_a],[mu_b,v_scales(i)*v_b],'bPlot',false);
    bd_shift=results_pretend.norm_bd;
    dprimes(i)=results_pretend.norm_dprime;
    results_shift=classify_normals(absent,present,'type','samp','reg',bd_shift,'samp_opt',false,'bPlot',false);
    
    samp_p11(i)=results_shift.samp_errmat(1,1)/n_samp;
    samp_p22(i)=results_shift.samp_errmat(2,2)/n_samp;
    
    norm_p11(i)=results_shift.norm_errmat(1,1)/.5;
    norm_p22(i)=results_shift.norm_errmat(2,2)/.5;    
end
norm_p11_sd=sqrt(norm_p11.*(1-norm_p11)/n_samp);
norm_p22_sd=sqrt(norm_p22.*(1-norm_p22)/n_samp);


figure; hold on
xline(dprime,'k',"actual d'",'LabelVerticalAlignment','bottom') % optimal criterion

% p11
x=dprimes;
y=norm_p11;
dy=norm_p11_sd;
fill([x;flipud(x)],[y-dy;flipud(y+dy)],[0 0 1],'facealpha',.3,'linestyle','none');
plot(dprimes,samp_p11,'.b','markersize',6)

% p22
y=norm_p22;
dy=norm_p22_sd;
fill([x;flipud(x)],[y-dy;flipud(y+dy)],[1 0 0],'facealpha',.3,'linestyle','none');
plot(dprimes,samp_p22,'.r','markersize',6)

xlim([min(dprimes) 10]); ylim([0 1])
xlabel('assumed d''')
set(gca,'fontsize',13); box off