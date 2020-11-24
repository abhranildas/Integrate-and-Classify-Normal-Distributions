load camouflage_edge_data
absent=[edge_powers_2(:,1),edge_lpr_2(:,1),edge_powers_4(:,1),edge_lpr_4(:,1),edge_powers_8(:,1),edge_lpr_8(:,1)];
present=[edge_powers_2(:,2),edge_lpr_2(:,2),edge_powers_4(:,2),edge_lpr_4(:,2),edge_powers_8(:,2),edge_lpr_8(:,2)];
N_samp=size(absent,1);
n_samp=1e2;
n_bs=1e3; % # of bootstrap samples

results=classify_normals(absent,present,'type','samp','samp_opt',false,'bPlot',false);
bd=results.norm_bd;

q0s=linspace(250,410,100)';

true_p1_1=nan(size(q0s));
true_p2_2=nan(size(q0s));

norm_p1_1=nan(size(q0s));
norm_p2_2=nan(size(q0s));

norm_l_mu=nan(size(q0s));
norm_l_v=nan(size(q0s));

true_bs_p1_1=nan(length(q0s),n_bs);
true_bs_p2_2=nan(length(q0s),n_bs);

true_bs_l=nan(length(q0s),n_bs);

parfor i=1:length(q0s)
    i
    bd_shift=bd;
    bd_shift.q0=q0s(i);
    results_shift=classify_normals(absent,present,'type','samp','dom',bd_shift,'samp_opt',false,'bPlot',false);
    
    % true relative fractions
    true_p1_1(i)=results_shift.samp_errmat(1,1)/N_samp;
    true_p2_2(i)=results_shift.samp_errmat(2,2)/N_samp;
    
    % normal outcome probabilities
    norm_p1_1(i)=results_shift.norm_errmat(1,1)/.5;
    norm_p2_2(i)=results_shift.norm_errmat(2,2)/.5;
    
    % normal likelihood sampling mean and sd
    norm_l11_mu=sum(binopdf(0:n_samp,n_samp,norm_p1_1(i)).^2);
    norm_l22_mu=sum(binopdf(0:n_samp,n_samp,norm_p2_2(i)).^2);
    
%     norm_l11_v=sum(binopdf(0:n_samp,n_samp,norm_p1_1(i)).^3)-norm_l11_mu^2;
%     norm_l22_v=sum(binopdf(0:n_samp,n_samp,norm_p2_2(i)).^3)-norm_l22_mu^2;
    
    norm_l_mu(i)=norm_l11_mu*norm_l22_mu;    
    norm_l_v(i)=sum(binopdf(0:n_samp,n_samp,norm_p1_1(i)).^3)...
               *sum(binopdf(0:n_samp,n_samp,norm_p2_2(i)).^3)...
               -norm_l11_mu^2*norm_l22_mu^2;
%     norm_l_v(i)=(norm_l11_v+norm_l11_mu^2)*(norm_l22_v+norm_l22_mu^2)-norm_l11_mu^2*norm_l22_mu^2;
    
    for j=1:n_bs
        % classify samples from the true distributions
        results_shift_bs=classify_normals(absent((j-1)*n_samp+1:j*n_samp,:),...
            present((j-1)*n_samp+1:j*n_samp,:),'type','samp','dom',bd_shift,...
            'samp_opt',false,'bPlot',false);
        
        % outcome counts
        samp_c11=results_shift_bs.samp_errmat(1,1);
        samp_c22=results_shift_bs.samp_errmat(2,2);
        
        % relative fractions
        true_bs_p1_1(i,j)=samp_c11/n_samp;
        true_bs_p2_2(i,j)=samp_c22/n_samp;
        
        % error matrix likelihood under normal
        true_bs_l(i,j)=binopdf(samp_c11,n_samp,norm_p1_1(i))*...
            binopdf(samp_c22,n_samp,norm_p2_2(i));        
        
%         % classify samples from the normal distributions
%         results_shift_bs=classify_normals(mvnrnd(mean(absent),cov(absent),n_samp),...
%             mvnrnd(mean(present),cov(present),n_samp),'type','samp','dom',bd_shift,...
%             'samp_opt',false,'bPlot',false);
%         
%         % outcome counts
%         samp_c11=results_shift_bs.samp_errmat(1,1);
%         samp_c22=results_shift_bs.samp_errmat(2,2);
%         
%         % error matrix likelihood under normal
%         norm_bs_l(i,j)=binopdf(samp_c11,n_samp,norm_p1_1(i))*...
%             binopdf(samp_c22,n_samp,norm_p2_2(i));
    end
end
norm_p1_1_sd=sqrt(norm_p1_1.*(1-norm_p1_1)/n_samp);
norm_p2_2_sd=sqrt(norm_p2_2.*(1-norm_p2_2)/n_samp);

figure; hold on
xline(bd.q0,'k','optimal boundary','LabelVerticalAlignment','bottom') % optimal criterion

% p11
x=q0s;
y=true_p1_1;
dy=std(true_bs_p1_1,0,2);
fill([x;flipud(x)],[y-dy;flipud(y+dy)],[0.9290 0.6940 0.1250],'facealpha',.3,'linestyle','none');

y=norm_p1_1;
dy=norm_p1_1_sd;
fill([x;flipud(x)],[y-dy;flipud(y+dy)],'k','facecolor','none','edgecolor',[0.9290 0.6940 0.1250]);

% p22
y=true_p2_2;
dy=std(true_bs_p2_2,0,2);
fill([x;flipud(x)],[y-dy;flipud(y+dy)],[0.6 0.4 0.2],'facealpha',.3,'linestyle','none');

y=norm_p2_2;
dy=norm_p2_2_sd;
fill([x;flipud(x)],[y-dy;flipud(y+dy)],'k','facecolor','none','edgecolor',[0.6 0.4 0.2]);

xlim([min(q0s) max(q0s)]); ylim([-.01 1.01])
xlabel('boundary offset q_0')
set(gca,'fontsize',13); box off

% sample likelihood
figure; hold on
xline(bd.q0,'k','optimal boundary','LabelVerticalAlignment','bottom') % optimal criterion
y=mean(true_bs_l,2);
dy=std(true_bs_l,0,2);
fill([x;flipud(x)],[y-dy;flipud(y+dy)],'k','facealpha',.3,'linestyle','none');

y=norm_l_mu;
dy=sqrt(norm_l_v);
fill([x;flipud(x)],[y-dy;flipud(y+dy)],'k','facecolor','none','edgecolor','k');

xlim([min(q0s) max(q0s)]); ylim([-.1 1.15])
xlabel('boundary offset q_0')
set(gca,'fontsize',13); box off

%%
% normals=struct;
% normals(1).mu=[1;0]; normals(1).v=2*eye(2);
% normals(2).mu=[0;1]; normals(2).v=eye(2);
% normals(3).mu=[-1;0]; normals(3).v=eye(2);
% normals(4).mu=[0;-1]; normals(4).v=eye(2);


N_samp=1e4;
n_samp=1e2;
n_bs=100; % # of bootstrap samples

normals=struct;
normals(1).mu=[0;-4;-2;-2]; normals(1).v=diag([1 2 3 4]);
normals(2).mu=[1;1;0;0]; normals(2).v=eye(4)/4;
normals(3).mu=[3;4;5;7]; normals(3).v=eye(4);
normals(4).mu=[5;5;7;4]; normals(4).v=eye(4);

priors=[1 1 1 1]/4;

samples=struct;
for i=1:4
    samples(i).sample=mvtrnd(normals(i).v,3,N_samp)+normals(i).mu';
end

% results_samp=classify_normals_multi(samples,'type','samp','priors',priors,'mc_samples',1e3,'plotmode',[1;1;1;1])

% now classify using samples from t distributions similar to these normals
% samples=struct;
% for i=1:4
%     samples(i).sample=mvtrnd(normals(i).v,3,1e4)+normals(i).mu';
% end

prior_fam=10.^linspace(-6,8,30)';

true_p1_1=nan(size(prior_fam));
true_pe=nan(size(prior_fam));

norm_p1_1=nan(size(prior_fam));
norm_pe=nan(size(prior_fam));
norm_pe_v=nan(size(prior_fam));
norm_bs_pe=nan(length(prior_fam),n_bs);

true_bs_p1_1=nan(length(prior_fam),n_bs);
true_bs_pe=nan(length(prior_fam),n_bs);

parfor i=1:length(prior_fam)
    priors_shift=[1 prior_fam(i) 1 1];
    priors_shift=priors_shift/sum(priors_shift);
    domains=cell(length(normals),1);
    for j=1:length(normals)
        domains{j}=@(n,mu,v) opt_class_multi(n,normals,j,'mu',mu,'v',v,'priors',priors_shift);
    end
    
    results_shift=classify_normals_multi(samples,'type','samp','doms',domains,'plotmode',false);
    
    % true outcomes
    true_p1_1(i)=results_shift.samp_errmat(2,2)/N_samp;
    true_pe(i)=results_shift.samp_err;
    
    % normal outcomes
    norm_p1_1(i)=results_shift.norm_errmat(2,2)/sum(results_shift.norm_errmat(2,:));
    norm_pe(i)=results_shift.norm_err;
    errs=sum(~eye(length(samples)).*results_shift.norm_errmat,2);
%     errs=errs./priors';
%     norm_pe_v(i)=sum(priors'.^2.*errs.*(1-errs))/n_samp;
    norm_pe_v(i)=sum(errs.*(priors'-errs))/n_samp;
    
    for j=1:n_bs
        [i j]
        % classify samples from the true distributions
        samples_bs=struct;
        for k=1:length(samples)
            samples_bs(k).sample=datasample(samples(k).sample,n_samp);
        end
        results_shift_bs=classify_normals_multi(samples_bs,'type','samp','doms',domains,'plotmode',false);
        
        % error matrix
        true_bs_p1_1(i,j)=results_shift_bs.samp_errmat(2,2)/n_samp;
        % overall error
        true_bs_pe(i,j)=results_shift_bs.samp_err;        
        
        % classify samples from the normal distributions
%         samples_bs=struct;
%         for k=1:length(samples)
%             samples_bs(k).sample=mvnrnd(normals(k).mu,normals(k).v,n_samp);
%         end
% 
%         results_shift_bs=classify_normals_multi(samples_bs,'type','samp','doms',domains,'plotmode',false);
%         norm_bs_pe(i,j)=results_shift_bs.samp_err;
        
    end
end
norm_p1_1_sd=sqrt(norm_p1_1.*(1-norm_p1_1)/n_samp);

figure; hold on
xline(priors(1),'k','optimal boundary','LabelVerticalAlignment','bottom') % optimal criterion

% p11
x=prior_fam/3;%./(prior_fam+sum(priors(2:end)));
y=true_p1_1;
dy=std(true_bs_p1_1,0,2);
fill([x;flipud(x)],[y-dy;flipud(y+dy)],[0.9290 0.6940 0.1250],'facealpha',.3,'linestyle','none');

y=norm_p1_1;
dy=norm_p1_1_sd;
fill([x;flipud(x)],[y-dy;flipud(y+dy)],'k','facecolor','none','edgecolor',[0.9290 0.6940 0.1250]);

% pe
y=true_pe;
dy=std(true_bs_pe,0,2);
fill([x;flipud(x)],[y-dy;flipud(y+dy)],[0.6 0.4 0.2],'facealpha',.3,'linestyle','none');

y=norm_pe;
dy=sqrt(norm_pe_v);
% dy=std(norm_bs_pe,0,2);
fill([x;flipud(x)],[y-dy;flipud(y+dy)],'k','facecolor','none','edgecolor',[0.6 0.4 0.2]);

xlim([min(x) max(x)]); set(gca,'xscale','log'); ylim([0 1])
xlabel('assumed prior p_1')
set(gca,'fontsize',13); box off

