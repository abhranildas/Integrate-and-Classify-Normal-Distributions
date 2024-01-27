function p_ray=prob_ray(init_sign_ray,z_ray,dim,varargin)

    % parse inputs
    parser=inputParser;
    parser.KeepUnmatched=true;
    addRequired(parser,'init_sign_ray',@isnumeric);
    addRequired(parser,'z_ray',@isnumeric);
    addRequired(parser,'dim',@isnumeric);
    addParameter(parser,'side','normal',@(x) strcmpi(x,'normal') || strcmpi(x,'complement'));
    addParameter(parser,'vpa',false,@islogical);

    parse(parser,init_sign_ray,z_ray,dim,varargin{:});
    side=parser.Results.side;
    vpaflag=parser.Results.vpa;

    % function to compute probability slice along each ray
    [f_big,f_small]=Phibar_ray_split(z_ray,dim);
    p_ray_big=init_sign_ray+1+init_sign_ray*sum((-1).^(1:length(z_ray)).*f_big);
    p_ray_small=init_sign_ray*sum((-1).^(1:length(z_ray)).*f_small);
    if strcmpi(side,'normal')
        p_ray=p_ray_big+p_ray_small;
    elseif strcmpi(side,'complement')
        p_ray=2-p_ray_big-p_ray_small;
    end

    % if there are roots but numeric answer is 0, evaluate symbolically
    if vpaflag&&numel(z_ray)&&(~p_ray)
        p_ray=prob_ray_sym(init_sign_ray,z_ray,dim,side);
    end

end