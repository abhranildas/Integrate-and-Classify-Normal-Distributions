function [p,bd_pts]=int_norm_ray(mu,v,dom,varargin)
    % Integrate a multinormal distribution over a specified domain, using
    % the ray method.

    parser = inputParser;
    parser.KeepUnmatched=true;
    addRequired(parser,'mu',@isnumeric);
    addRequired(parser,'v',@isnumeric);
    addRequired(parser,'dom');
    addParameter(parser,'dom_type','quad');
    addParameter(parser,'output','prob'); % probability or probability density
    addParameter(parser,'force_mc',false,@islogical);
    addParameter(parser,'fun_span',5);
    addParameter(parser,'fun_resol',100);
    addParameter(parser,'AbsTol',1e-10);
    addParameter(parser,'RelTol',1e-2);
    addParameter(parser,'mc_samples',500);

    parse(parser,mu,v,dom,varargin{:});
    output=parser.Results.output;
    force_mc=parser.Results.force_mc;
    AbsTol=parser.Results.AbsTol;
    RelTol=parser.Results.RelTol;
    mc_samples=parser.Results.mc_samples;

    dim=length(mu);

    global bd_pts
    bd_pts=[];

    if force_mc || dim>3
        % Monte-Carlo integration

        % uniform random rays (points on n-sphere)
        n_z=mvnrnd(zeros(dim,1),eye(dim),mc_samples)';
        n_z=n_z./vecnorm(n_z,2);

        p_rays=norm_prob_across_rays(mu,v,dom,n_z,varargin{:});

        if iscell(p_rays) % if p_rays is a cell array
            % separate numeric and symbolic
            num_idx=cellfun(@isnumeric, p_rays);
            sum_num=sum(horzcat(p_rays{num_idx}));
            sum_sym=sum(horzcat(p_rays{~num_idx}));

            if logical(vpa(sum_sym)~=0) && (double(log10(sum_sym)) > log10(sum_num) + log10(RelTol))
                p=(sum_num+sum_sym)/numel(p_rays);
                warning('Probability too small for double precision. Returning as symbol, use vpa to evaluate.')
            else
                if logical(vpa(sum_sym)==0)
                    warning('Probability on some rays too small even for variable precision. Returning 0.')
                end
                n_num=sum(num_idx); % number of numeric rays
                if n_num==0
                    p=0;
                else
                    p=sum_num/sum(num_idx);
                end
                %             [log10(sum_num) double(log10(sum_sym))]
                % incorporate the symbolic only if it's
                % bigger than RelTol of the numeric
            end
        else
            p=mean(p_rays);
        end

        if strcmpi(output,'prob')
            p=p/2;
        end

    else
        % grid integration

        if dim==1
            p=norm_prob_across_angles(mu,v,dom,varargin{:});
        elseif dim==2
            % p=integral(@(theta) norm_prob_across_angles(mu,v,dom,'theta',theta,varargin{:}),0,pi,'AbsTol',AbsTol,'RelTol',RelTol);
            syms theta
            p=vpaintegral(@(theta) norm_prob_across_angles(mu,v,dom,'theta',theta,varargin{:}),theta,0,pi,'AbsTol',AbsTol,'RelTol',RelTol);
        elseif dim==3
            p=integral2(@(theta,phi) norm_prob_across_angles(mu,v,dom,'theta',theta,'phi',phi,varargin{:}),0,pi/2,0,2*pi,'AbsTol',AbsTol,'RelTol',RelTol);
        end
    end