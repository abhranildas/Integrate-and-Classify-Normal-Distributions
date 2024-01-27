function p=int_norm_quad_gx2(mu,v,quad,varargin)
    % Find the probability that a quadratic form of a normal variate x
    % x'q2x + q1'x + q0 >= 0
    % using the generalized chi-squared CDF (Imhof's method).
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

    parser = inputParser;
    parser.KeepUnmatched=true;
    addRequired(parser,'mu',@isnumeric);
    addRequired(parser,'v',@isnumeric);
    addRequired(parser,'quad');
    addParameter(parser,'AbsTol',1e-10);
    addParameter(parser,'RelTol',1e-2);
    addParameter(parser,'side','normal',@(x) strcmpi(x,'normal') || strcmpi(x,'complement') );
    parse(parser,mu,v,quad,varargin{:});

    side=parser.Results.side;

    if ~nnz(quad.q2) % if q2 is zero, linear discriminant
        quad_s=standard_quad(quad,mu,v); % standardize quad
        if strcmpi(side,'normal')
            p=normcdf(quad_s.q0/norm(quad_s.q1));
        elseif strcmpi(side,'complement')
            p=normcdf(quad_s.q0/norm(quad_s.q1),'upper');
        end
        %     pc=normcdf(-quad_s.q0/norm(quad_s.q1));
    else
        % get generalized chi-squared parameters
        [w,k,lambda,m,s]=norm_quad_to_gx2_params(mu,v,quad);
        %     pc=gx2cdf(0,w,k,lambda,m,s,varargin{:});
        if strcmpi(side,'normal')
            p=gx2cdf(0,w,k,lambda,m,s,'upper',varargin{:});
        elseif strcmpi(side,'complement')
            p=gx2cdf(0,w,k,lambda,m,s,varargin{:});
        end
        %     p=gx2cdf(0,w,k,lambda,m,s,'upper',varargin{:});
    end
