mu_2d_1=[0;0;0];
v_3d_1=eye(3);

mu_3d_2=[1;1;1];
v_3d_2=[4    0  -.3  ;
    0    1  -.5  ;
    -.3  -.5    2  ];

sigma=10.^linspace(-6,2,100);
alist=10.^linspace(.01,3,7);

d_1d_ovr=nan(length(alist),length(sigma));
d_1d_rms=nan(size(d_1d_ovr));
d_1d_avg=nan(size(d_1d_ovr));

d_3d_ovr=nan(size(d_1d_ovr));
d_3d_rms=nan(size(d_1d_ovr));
d_3d_avg=nan(size(d_1d_ovr));

for i=1:length(alist)
    a=alist(i);
    parfor j=1:length(sigma)
        [i j]
        results_1d=classify_normals([0,(a*sigma(j))^2],[1,sigma(j)^2],'plotmode',false);
        d_1d_ovr(i,j)=results_1d.norm_dprime;
        
        results_3d=classify_normals([mu_2d_1,(a*sigma(j))^2*v_3d_1],[mu_3d_2,sigma(j)^2*v_3d_2],'plotmode',false);
        d_3d_ovr(i,j)=results_3d.norm_dprime;
        
        v_3d_rms=((a*sigma(j))^2*v_3d_1+sigma(j)^2*v_3d_2)/2;
        d_3d_rms(i,j)=sqrt((mu_2d_1-mu_3d_2)'/v_3d_rms*(mu_2d_1-mu_3d_2));
        
        v_3d_avg=((a*sigma(j)*sqrtm(v_3d_1)+sigma(j)*sqrtm(v_3d_2))/2)^2;
        d_3d_avg(i,j)=sqrt((mu_2d_1-mu_3d_2)'/v_3d_avg*(mu_2d_1-mu_3d_2));
    end
    d_1d_rms(i,:)=1./(sigma*sqrt((a^2+1)/2));
    d_1d_avg(i,:)=2./(sigma*(a+1));
end

figure;
subplot(1,2,1); hold on; axis([0 70 0 1.1])
ylabel("$d'_a/d'_o, d'_e/d'_o$",'interpreter','latex')
set(gca,'fontsize',13,'xtick',[0 70])
subplot(1,2,2); hold on; axis([0 70 0 1.1])
set(gca,'fontsize',13,'xtick',[0 70])
cmap=colormap(flipud(winter(7))); c=colorbar('ticks',linspace(0,1,4),'ticklabels',[1 10 100 1000]);

for i=1:size(d_1d_ovr,1)
    subplot(1,2,1)
    plot(d_1d_ovr(i,:),d_1d_avg(i,:)./d_1d_ovr(i,:),'color',cmap(i,:))
    plot(d_1d_ovr(i,:),d_1d_rms(i,:)./d_1d_ovr(i,:),'color',cmap(i,:))
    
    subplot(1,2,2)
    plot(d_3d_ovr(i,:),d_3d_avg(i,:)./d_3d_ovr(i,:),'color',cmap(i,:))
    plot(d_3d_ovr(i,:),d_3d_rms(i,:)./d_3d_ovr(i,:),'color',cmap(i,:))
end

%%
mu_a=0; sigma=.05; a=2;
mu_b=1; s_b=.05;
results=classify_normals([mu_a,(a*sigma)^2],[mu_b,sigma^2]);
axis([-.5 2.5 0 4.2])
title ''
set(gca,'xtick',[0 1],'ytick',[],'box','off','fontsize',13,'tickdir','out')
%%
mu_n=0; mu_s=2;
v_n=0.5; v_s=1;

mu_s_list=linspace(0,2,100);

d_o=nan(size(mu_s_list));
d_a=nan(size(mu_s_list));
d_e=nan(size(mu_s_list));

for i=1:length(mu_s_list)
    mu_s=mu_s_list(i)
    
    mu_a=[mu_s;mu_n];
    v_a=[v_s 0; 0 v_n];
    
    mu_b=[mu_n;mu_s];
    v_b=[v_n 0; 0 v_s];
    
    results=classify_normals([mu_a,v_a],[mu_b,v_b],'method','gx2','AbsTol',0,'RelTol',1e-5,'plotmode',false);
    % axis image; axis(2*mu_s*[-1 1 -1 1])
    d_o(i)=results.norm_dprime;
    
%     v_rms=(v_a+v_b)/2;
%     results=classify_normals([mu_a,v_rms],[mu_b,v_rms],'plotmode',false);
    % axis image; axis(2*mu_s*[-1 1 -1 1])
    d_a(i)=results.norm_dprime_a;
    
%     v_avg=((sqrtm(v_a)+sqrtm(v_b))/2)^2;
%     results=classify_normals([mu_a,v_avg],[mu_b,v_avg],'plotmode',false);
    % axis image; axis(2*mu_s*[-1 1 -1 1])
    d_e(i)=results.norm_dprime_e;
end

figure; hold on
plot(mu_s_list,d_a./d_o);
plot(mu_s_list,d_e./d_o);