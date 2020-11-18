function samp_errmat=samp_errmat_multi(normals,samples,varargin)
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
n_dists=length(normals);
parser=inputParser;
parser.KeepUnmatched=true;
addRequired(parser,'normals',@isstruct);
addRequired(parser,'samples',@isstruct);
addParameter(parser,'priors',ones(1,n_dists)/n_dists,@(x) isnumeric(x) && all(x > 0) && all(x < 1));
addParameter(parser,'vals',eye(n_dists), @(x) isnumeric(x) && ismatrix(x));

parse(parser,normals,samples,varargin{:});
priors=parser.Results.priors;
vals=parser.Results.vals;
val_gains=dot((eye(n_dists)-ones(n_dists)/2)*2,vals,2);
samp_errmat=nan(n_dists);
for i_samp=1:n_dists
    log_ex_val_gains=nan(size(samples(i_samp).sample,1),n_dists);
    for i_norm=1:n_dists
        log_ex_val_gains(:,i_norm)=log_ex_val_gain(samples(i_samp).sample,normals(i_norm).mu,normals(i_norm).v,priors(i_norm),val_gains(i_norm));
    end
    [~,class_selected]=max(log_ex_val_gains,[],2);
    samp_errmat(i_samp,:)=histcounts(class_selected,1:n_dists+1);
end
end


