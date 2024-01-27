function [p_rays,bd_pts_rays]=norm_prob_across_rays(mu,v,dom,n_z,varargin)

    % parse inputs
    parser=inputParser;
    parser.KeepUnmatched=true;
    addRequired(parser,'mu',@isnumeric);
    addRequired(parser,'v',@isnumeric);
    addRequired(parser,'dom',@(x) isstruct(x) || isa(x,'function_handle') || ismatrix(x));
    addParameter(parser,'output','prob'); % probability or probability density
    addParameter(parser,'side','normal');
    addParameter(parser,'dom_type','quad');
    addParameter(parser,'fun_span',5);
    addParameter(parser,'fun_resol',100);
    addParameter(parser,'fun_grad',[],@(x) isa(x,'function_handle'));
    addParameter(parser,'n_bd_pts',1e4);
    addParameter(parser,'add_bd_pts',true);

    parse(parser,mu,v,dom,varargin{:});

    output=parser.Results.output;
    dim=length(mu);

    % ray-trace the standardized domain/function
    dom_standard_raytrace=@(n) standard_ray_trace(dom,n,'mu',mu,'v',v,varargin{:});

    if strcmpi(output,'prob')

        % initial signs and boundary distances in standardized space
        [init_sign,z]=dom_standard_raytrace(n_z);

        % probability on rays
        p_rays=cellfun(@(init_sign_ray,z_ray) prob_ray(init_sign_ray,z_ray,dim,varargin{:}),num2cell(init_sign),z,'un',0);

    elseif strcmpi(output,'prob_dens') % probability density calculations

        % roots of standardized function along rays
        [~,z]=dom_standard_raytrace(n_z);

        % gradient of standardized function
        fun_grad=parser.Results.fun_grad;
        standard_gradf=@(z) standard_fun_grad(z,fun_grad,mu,v);

        % probability density on rays
        p_rays=cellfun(@(n_ray,z_ray) prob_dens_ray(n_ray,z_ray,standard_gradf), num2cell(n_z,1),z);
    end

    % if p_rays is a cell but there are no symbols, convert to numeric array
    if iscell(p_rays)
        num_idx=cellfun(@isnumeric, p_rays);
        if all(num_idx)
            p_rays=cell2mat(p_rays);
        end
    end

    % standard boundary points
    std_bd_pts_ray=cellfun(@(z_ray,n_ray) z_ray.*n_ray, z,num2cell(n_z,1),'un',0);

    % boundary points
    bd_pts_rays=sqrtm(v)*horzcat(std_bd_pts_ray{:})+mu;
    if parser.Results.add_bd_pts
        global bd_pts
        bd_pts=[bd_pts,bd_pts_rays];
    end


