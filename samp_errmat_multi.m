function [samp_errmat,samp_class]=samp_errmat_multi(samples,doms)
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

% parse inputs
parser=inputParser;
parser.KeepUnmatched=true;
addRequired(parser,'samples',@isstruct);
addRequired(parser,'doms',@iscell);
parse(parser,samples,doms);

n_dists=length(samples);
samp_class=cell(n_dists,1);
samp_errmat=nan(n_dists);
for i_samp=1:n_dists
    samp=samples(i_samp).sample;
    samp_class_this=nan(size(samp,1),1);
    for i_norm=1:n_dists
        dom=doms{i_norm};
        [~,~,samp_correct_this]=dom(samp',[]);
        samp_class_this(samp_correct_this)=i_norm;
        samp_errmat(i_samp,i_norm)=nnz(samp_correct_this);
    end
    samp_class{i_samp}=samp_class_this;
end
end


