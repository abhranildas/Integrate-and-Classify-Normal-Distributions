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
    elseif isvector(plotmode)
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

norm_bd_pts=cell(n_dists,1);
norm_errmat=nan(n_dists);
norm_errs=nan(n_dists,1);
% compute accuracy and boundary for each normal, and combine
for i=1:n_dists % integrate each normal
    for j=1:n_dists % within the boundary of each normal
        [p,pc,bd_pts]=integrate_normal(normals(i).mu,normals(i).v,...
            @(n,mu,v) opt_class_multi(n,normals,j,'mu',mu,'v',v,varargin{:}),...
            'reg_type','ray_scan','bPlot',false,varargin{:});
        fprintf('Integrating normal %d in region %d\n',[i j])
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
    samp_counts=samp_counts_multi(normals,dists,varargin{:});
    results.samp_errmat=samp_counts;
    %results.samp_errs=samp_counts./sum(samp_counts,2);
    samp_err=sum(samp_counts(~eye(n_dists)))/sum(samp_counts(:));
    results.samp_err=samp_err;
    if ~isequal(vals,eye(n_dists)) % if outcome values are supplied
        results.samp_valmat=samp_counts.*vals;
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

    
    % Colormap of the different regions
    %     fcontour(@(x,y) multifun(x,y,normals,1),'levellist',0.5:1:4.5,'fill','on','linecolor','none');
    %     caxis([.5 4.5]);
    %     colormap((colors(1:4,:)+[1 1 1])/2)
    
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
            plot_boundary(@(n,mu,v) opt_class_multi(n,normals,i,'mu',mu,'v',v,varargin{:}),dim,'mu',normals(i).mu,'reg_type','ray_scan')
        end
    end
    hold off
end