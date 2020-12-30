function results=classify_normals_multi(dists,varargin)
	% CLASSIFY_NORMALS_MULTI Compute quantities concerning classification
	% accuracy between multiple normal distributions, or samples.
	%
	% Abhranil Das <abhranil.das@utexas.edu>
	% Center for Perceptual Systems, University of Texas at Austin
	% If you use this code, please cite:
	% <a href="matlab:web('https://arxiv.org/abs/2012.14331')"
	% >A method to integrate and classify normal distributions</a>.
	%
    % Example:
    % normals=struct;
    % normals(1).mu=[1;0]; normals(1).v=2*eye(2);
    % normals(2).mu=[0;1]; normals(2).v=eye(2);
    % normals(3).mu=[-1;0]; normals(3).v=eye(2);
    % normals(4).mu=[0;-1]; normals(4).v=eye(2);
    % results=classify_normals_multi(normals)
    %
	% Required inputs:
	% dists         struct of parameters, with each element containing mu and v
	%               or struct of samples, with each element containing a sample
	%
	% Optional name-value inputs:
	% doms          custom (non-optimal) classification domains, as a cell
	%               array, whose elements are ray-scan functions defining the
	%               the classification domain of each normal
	% priors        vector of priors of the normals. Default is equal priors.
	% vals          matrix of outcome values. v(i,j) is the value of
	%               classifying a sample from i as j.
	% input_type    'norm' for normal parameter inputs (default), 'samp' for
	%               sample inputs.
	% AbsTol        absolute tolerance for the error rate computations.
	% RelTol        relative tolerance for the error rate computations.
	%               The absolute OR the relative tolerance will be satisfied.
	% plotmode      true (default) for plot, or false, or column vector along
	%               which a >3d multi-class problem will be projected (default
	%               is along the first dimension.
	%
	% Output: struct containing
	% norm_bd_pts   cell array containing points on the optimal boundaries
	%               around each normal computed by the ray-scan integration
	% norm_errmat   error matrix. e(i,j)=prob. of classifying a sample from
	%               normal i as j.
	% norm_err      overall error rate
	% norm_valmat   matrix of expected classification outcome values. v(i,j)=
	%               expected value of classifying sample from i as j. Returned
	%               only when custom values are supplied.
	% norm_val      overall expected outcome value. Returned only when custom
	%               values are suppliued.
	% samp_errmat   error matrix of classified samples. e(i,j)= no. of i
	%               samples classified as j.
	% samp_err      overall error rate of classifying the samples.
	% samp_valmat   matrix of total sample classification outcome values.
	%               Returned only when custom values are supplied.
	% samp_val      total value of classifying the samples. Returned only when
	%               custom values are supplied.
	%
	% See also:
	% <a href="matlab:open(strcat(fileparts(which('integrate_normal')),filesep,'doc',filesep,'GettingStarted.mlx'))">Interactive demos</a>
	% classify_normals
	
	% parse inputs
	n_dists=length(dists);
	parser=inputParser;
	parser.KeepUnmatched=true;
	addRequired(parser,'dists',@isstruct);
	addParameter(parser,'doms',[],@iscell);
	addParameter(parser,'priors',ones(1,n_dists)/n_dists,@(x) isnumeric(x) && all(x > 0) && all(x < 1));
	addParameter(parser,'vals',eye(n_dists),@(x) isnumeric(x) && ismatrix(x));
	addParameter(parser,'input_type','norm',@(s) strcmpi(s,'norm') || strcmpi(s,'samp'));
	addParameter(parser,'plotmode',true, @(x) islogical(x) || iscolumn(x));
	
	parse(parser,dists,varargin{:});
	if strcmpi(parser.Results.input_type,'samp')
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
	if strcmpi(parser.Results.input_type,'samp') && any(strcmpi(parser.UsingDefaults,'priors'))
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
				doms{j},'dom_type','ray_scan',varargin{:},'plotmode',false,varargin{:});
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
	if strcmpi(parser.Results.input_type,'samp')
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
		if strcmpi(parser.Results.input_type,'norm')
			title(sprintf("error = %g",norm_err),'interpreter','latex')
		elseif strcmpi(parser.Results.input_type,'samp')
			title(sprintf("error = %g / %g",[norm_err,samp_err]))
		end
		
		% plot samples
		if strcmpi(parser.Results.input_type,'samp')
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