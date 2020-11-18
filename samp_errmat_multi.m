function samp_errmat=samp_errmat_multi(samples,doms)
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
samp_errmat=nan(n_dists);
for i_samp=1:n_dists
    for i_norm=1:n_dists
        dom=doms{i_norm};
        [~,~,samp_correct]=dom(samples(i_samp).sample',[],[]);
        samp_errmat(i_samp,i_norm)=nnz(samp_correct);
    end
end
end


