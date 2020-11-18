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
parser=inputParser;
parser.KeepUnmatched=true;
addParameter(parser,'doms',[],@iscell);
addRequired(parser,'dists',@isstruct);
addParameter(parser,'priors',ones(1,n_dists)/n_dists,@(x) isnumeric(x) && all(x > 0) && all(x < 1));
addParameter(parser,'vals',eye(n_dists),@(x) isnumeric(x) && ismatrix(x));
addParameter(parser,'type','norm',@(s) strcmpi(s,'norm') || strcmpi(s,'samp'));
addParameter(parser,'plotmode',true, @(x) islogical(x) || iscolumn(x));

parse(parser,dists,varargin{:});
if strcmpi(parser.Results.type,'samp')
    normals=struct;
    for i=1:n_dists
        normals(i).mu=mean(dists(i).sample)';
        normals(i).v=cov(dists(i).sample);
    end
else
    normals=dists;
end
dim=length(normals(1).mu);

vals=parser.Results.vals;
plotmode=parser.Results.plotmode;
if dim>3
    if isequal(plotmode,true)
        plotmode=[1;zeros(dim-1,1)];
    elseif ~isequal(plotmode,false)
        plotmode=plotmode/norm(plotmode);
    end
end

% if sample input and priors are not specified,
if strcmpi(parser.Results.type,'samp') && any(strcmpi(parser.UsingDefaults,'priors'))
    % set priors according to input sample sizes
    priors=arrayfun(@(x) length(x.sample),dists);
    priors=priors/sum(priors);
else
    priors=parser.Results.priors;
end

% classification domains
doms=parser.Results.doms;
if isempty(doms) % no domains supplied
    % compute optimal domains
    doms=cell(length(normals),1);
    for i=1:length(normals)
        doms{i}=@(n,mu,v) opt_class_multi(n,normals,i,'mu',mu,'v',v,varargin{:});
    end
end

norm_bd_pts=cell(n_dists,1);
norm_errmat=nan(n_dists);
norm_errs=nan(n_dists,1);

% compute accuracy and boundary for each normal, and combine
for i=1:n_dists % integrate each normal
    for j=1:n_dists % within the boundary of each normal
        [p,pc,bd_pts]=integrate_normal(normals(i).mu,normals(i).v,...
            doms{j},'dom_type','ray_scan',varargin{:},'plotmode',false);
        fprintf('Integrating normal %d in domain %d\n',[i j])
        norm_errmat(i,j)=p*priors(i);
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
results.norm_errmat=norm_errmat;
results.norm_err=norm_err;

if ~isequal(vals,eye(n_dists)) % if outcome values are supplied
    results.norm_valmat=norm_errmat.*vals; % conditional expected values
    results.norm_val=sum(results.norm_valmat(:));
end

% sample classification errors
if strcmpi(parser.Results.type,'samp')
    samp_errmat=samp_errmat_multi(dists,doms);
    results.samp_errmat=samp_errmat;
    samp_err=sum(samp_errmat(~eye(n_dists)))/sum(samp_errmat(:));
    results.samp_err=samp_err;
    if ~isequal(vals,eye(n_dists)) % if outcome values are supplied
        results.samp_valmat=samp_errmat.*vals;
        results.samp_val=sum(results.samp_valmat(:));
    end
end

%% Plot data
if ~isequal(plotmode,false)
    figure; hold on
    colors=colororder;
    % plot title
    if strcmpi(parser.Results.type,'norm')
        title(sprintf("error = %g",norm_err),'interpreter','latex')
    elseif strcmpi(parser.Results.type,'samp')
        title(sprintf("error = %g / %g",[norm_err,samp_err]))
    end
    
    % plot samples
    if strcmpi(parser.Results.type,'samp')
        for i=1:n_dists
            if dim<=3
                plot_sample(dists(i).sample,priors(i),colors(i,:))
            else
                plot_sample(dists(i).sample*plotmode,priors(i),colors(i,:))
            end
        end
    end
    
    % plot normals and boundaries
    for i=1:n_dists
        if dim<=3
            plot_normal(normals(i).mu,normals(i).v,priors(i),colors(i,:))
        else
            plot_normal(plotmode'*normals(i).mu,plotmode'*normals(i).v*plotmode,priors(i),colors(i,:))
        end
        hold on
        if dim<=3
            plot_boundary(doms{i},dim,'mu',normals(i).mu,'dom_type','ray_scan')
        end
    end
    hold off
end