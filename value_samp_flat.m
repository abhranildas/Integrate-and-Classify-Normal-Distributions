function v=value_samp_flat(dim,x,obs_a,obs_b,vals)
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

bd_coeffs=struct;
bd_coeffs.a2=reshape(x(1:dim^2),[dim dim])';
bd_coeffs.a1=x(dim^2+1:dim^2+dim);
bd_coeffs.a0=x(end);

v=value_samp(obs_a,obs_b,bd_coeffs,vals);