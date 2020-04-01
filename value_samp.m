function [samp_ex_val,samp_val_mat]=value_samp(samp_1,samp_2,reg_coeffs,vals)
% Expected value given two samples and a boundary. If outcome values
% are not additionally specified, this is the classification accuracy.
% Credits:
%   Abhranil Das <abhranil.das@utexas.edu>
%	R Calen Walshe
%	Wilson S Geisler
%	Center for Perceptual Systems, University of Texas at Austin
% If you use this code, please cite:
%   A new method to compute classification error
%   https://jov.arvojournals.org/article.aspx?articleid=2750251

a2=reg_coeffs.a2;
a1=reg_coeffs.a1;
a0=reg_coeffs.a0;

%[~,~,x_sign_1]=ray_scan(reg_coeffs,'quad',samp_1',zeros(dim,1));
%[~,~,x_sign_2]=ray_scan(reg_coeffs,'quad',samp_1',zeros(dim,1));

samp_1_correct=dot(samp_1,samp_1*a2',2) + samp_1*a1 + a0 > 0;
samp_2_correct=dot(samp_2,samp_2*a2',2) + samp_2*a1 + a0 < 0;

samp_val_mat=zeros(2);
samp_val_mat(1,1)=sum(samp_1_correct)*vals(1,1);
samp_val_mat(1,2)=sum(~samp_1_correct)*vals(1,2);
samp_val_mat(2,2)=sum(samp_2_correct)*vals(2,2);
samp_val_mat(2,1)=sum(~samp_2_correct)*vals(2,1);

samp_ex_val=sum(samp_val_mat(:))/(size(samp_1,1)+size(samp_2,1));

end


