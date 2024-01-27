% mu_1=[2;4]; v_1=[1 1.5; 1.5 3];
% mu_2=[4;4]; v_2=[3 0; 0 1];
% 
% n_samp=1e3;
% samp_1=mvtrnd(v_1,1.5,n_samp)+mu_1';
% samp_2=mvtrnd(v_2,1.5,n_samp)+mu_2';

results=classify_normals([0 1],[1 1.01]);

absent=importdata('target_absent.txt',',',1);
present=importdata('target_present.txt',',',1);
samp_1=absent.data(:,[2 3]); samp_2=present.data(:,[2 3]);

results=classify_normals(samp_1,samp_2,'input_type','samp','samp_opt','step');

results_smooth=classify_normals(samp_1,samp_2,'input_type','samp','samp_opt','smooth');
