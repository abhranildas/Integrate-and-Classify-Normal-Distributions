function results=classify_normals_multi(dists,varargin)
% Classification between two distributions.
%
% How to use this command:
% See github readme at https://github.com/abhranildas/classify
%
% Credits:
%   Abhranil Das <abhranil.das@utexas.edu>
%	Wilson S Geisler
%	Center for Perceptual Systems, University of Texas at Austin
% If you use this code, please cite:
%   A new method to compute classification error
%   https://jov.arvojournals.org/article.aspx?articleid=2750251

% parse inputs
n_dists=length(dists);
parser = inputParser;
addRequired(parser,'dists',@isstruct);
addParameter(parser,'priors',ones(1,n_dists)/n_dists,@(x) isnumeric(x) && all(x > 0) && all(x < 1));
addParameter(parser,'vals',eye(n_dists),@(x) isnumeric(x) && ismatrix(x));
addParameter(parser,'type','norm',@(s) strcmp(s,'norm') || strcmp(s,'samp'));
addParameter(parser,'bPlot',true,@islogical);

parse(parser,dists,varargin{:});
vals=parser.Results.vals;
bPlot=parser.Results.bPlot;

if strcmp(parser.Results.type,'norm')
    dim=length(dists(1).mu);
    mus=nan(dim,n_dists);
    vs=nan(dim,dim,n_dists);
    for i=1:n_dists
        mus(:,i)=dists(i).mu;
        vs(:,:,i)=dists(i).v;
    end
else
    dim=size(dists(1).sample,2);
    mus=nan(dim,n_dists);
    vs=nan(dim,dim,n_dists);
    for i=1:n_dists
        mus(:,i)=mean(dists(i).sample);
        vs(:,:,i)=cov(dists(i).sample);
    end
end

if dim>3
    error('Multi-class classification can only handle 3 or fewer dimensions.')
end

% if sample input and priors are not specified,
if strcmp(parser.Results.type,'samp') && any(strcmp(parser.UsingDefaults,'priors'))
    % set priors according to input sample sizes
    priors=arrayfun(@(x) length(x.sample),dists);
    priors=priors/sum(priors);
else
    priors=parser.Results.priors;
end

if bPlot, figure; hold on; end

norm_bd_pts=cell(n_dists,1);
norm_err_mat=nan(n_dists);
norm_errs=nan(n_dists,1);
% compute accuracy and boundary for each normal, and combine
for i=1:n_dists % integrate each normal
    for j=1:n_dists % within the boundary of each normal
        [p,pc,bd_pts]=integrate_normal(mus(:,i),vs(:,:,i),...
            @(n,orig) opt_reg_multi(n,mus,vs,'idx',j,'priors',priors,'vals',vals,'orig',orig),...
            'reg_type','ray_scan','bPlot',false);
        fprintf('Integrating normal %d in region %d\n',[i j])        
        norm_err_mat(i,j)=p;
        if j==i
            norm_errs(i)=pc;
        end
        norm_bd_pts{j}=[norm_bd_pts{j},bd_pts];
    end
end

% trim to unique boundary points
for i=1:n_dists
    norm_bd_pts{i}=uniquetol(norm_bd_pts{i}',1e-12,'Byrows',true,'Datascale',1)';
end

norm_err=dot(priors,norm_errs);

results.norm_bd_pts=norm_bd_pts;
results.norm_err_mat=norm_err_mat;
results.norm_err=norm_err;

if ~isequal(vals,eye(n_dists)) % if outcome values are supplied
    results.norm_val_mat=results.norm_err_mat.*vals; % conditional expected values
    results.norm_val=sum(sum(results.norm_val_mat.*priors'));
end

%% Plot data
if bPlot
    hold on
    title(sprintf("error = %g",norm_err)) % plot title
    
    % plot samples
    colors=colororder;
    for i=1:n_dists
        if strcmp(parser.Results.type,'samp')
            plot_sample(dists(i).sample,[],priors(1),colors(i,:))
        end
    end
    % plot normals and boundaries
    for i=1:n_dists
        plot_normal(mus(:,i),vs(:,:,i),norm_bd_pts{i},priors(i),colors(i,:))
        hold on
    end
    hold off
end