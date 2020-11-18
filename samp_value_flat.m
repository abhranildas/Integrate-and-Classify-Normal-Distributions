function v=samp_value_flat(x,samp_1,samp_2,vals)
% Given observation samples and flattened boundary coefficients x, returns the expected value.
% This is maximized to yield the optimal boundary coefficients.
% Credits:
%   Abhranil Das <abhranil.das@utexas.edu>
%	R Calen Walshe
%	Wilson S Geisler
%	Center for Perceptual Systems, University of Texas at Austin
% If you use this code, please cite:
%   A new method to compute classification error
%   https://jov.arvojournals.org/article.aspx?articleid=2750251

dim=size(samp_1,2);
q2=zeros(dim);
q2(triu(true(dim)))=x(1:(dim^2+dim)/2);
q2=q2+triu(q2,1)';

quad=struct;
quad.q2=q2;
quad.q1=x(end-dim:end-1);
quad.q0=x(end);

v=samp_value(samp_1,samp_2,quad,'vals',vals);