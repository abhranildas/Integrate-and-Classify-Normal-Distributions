function f=norm_fun_pdf(x,mu,v,fun,varargin)
    % NORM_FUN_PDF Compute the pdf of a scalar function of a (multi)normal
    % variable.
    %
    % Abhranil Das <abhranil.das@utexas.edu>
    % Center for Perceptual Systems, University of Texas at Austin
    % If you use this code, please cite:
    % <a href="matlab:web('https://arxiv.org/abs/2012.14331')"
    % >A method to integrate and classify normal distributions</a>.
    %
    % Example:
    % mu=[1;2]; v=[1 0.5; 0.5 4];
    % fun=@(x,y) sin(x)+cos(y);
    % norm_fun_pdf(0.5,mu,v,fun)
    %
    % Required inputs:
    % x             point(s) to evaluate the cdf at
    % mu            normal mean as column vector
    % v             normal variance-covariance matrix
    % fun           scalar function of the normal, in one of two forms:
    %               • struct containing coefficients a2 (matrix), a1 (column
    %                 vector) and a0 (scalar) of a quadratic function:
    %                 x'*a2*x + a1'*x + a0 > 0
    %               • handle to a scalar-valued function
    %
    % Optional name-value inputs:
    % pdf_method    method to compute the pdf. 'diff' for numerically
    %               differentiating the cdf, or 'ray' for directly using
    %               the ray method (requires supplying fun_grad).
    % fun_span      trace radius (in Mahalanobis distance) for function.
    %               Default=5.
    % fun_resol     resolution of tracing (finding roots) of function.
    %               Default=100.
    % fun_grad      handle to a function that returns the scalar or vector 
    %               gradient of fun
    % dx            step-size for numerically differentiating cdf
    % AbsTol        absolute tolerance to compute the cdf in 'diff' method
    % RelTol        relative tolerance to compute the cdf in 'diff' method
    %               The absolute OR the relative tolerance will be satisfied.
    %
    % also supports some of the name-value inputs of integrate_normal.
    %
    % Outputs:
    % f             pdf
    %
    % See also:
    % <a href="matlab:open(strcat(fileparts(which('integrate_normal')),filesep,'doc',filesep,'GettingStarted.mlx'))">Interactive demos</a>
    % norm_fun_cdf
    % norm_fun_inv

    % parse inputs
    parser=inputParser;
    parser.KeepUnmatched=true;
    addRequired(parser,'x',@isnumeric);
    addRequired(parser,'mu',@isnumeric);
    addRequired(parser,'v',@isnumeric);
    addRequired(parser,'fun',@(x) isstruct(x)|| isa(x,'function_handle'));
    addParameter(parser,'pdf_method','ray');
    addParameter(parser,'dx',1e-3,@(x) isreal(x) && isscalar(x) && (x>=0));

    parse(parser,x,mu,v,fun,varargin{:});
    pdf_method=parser.Results.pdf_method;

    if isa(fun,'function_handle')
        if strcmpi(pdf_method,'ray')
            f=arrayfun(@(x_each) int_norm_ray(mu,v,fun,'dom_type','fun','output','prob_dens','fun_level',x_each,varargin{:}), x);
        elseif strcmpi(pdf_method,'diff')
            dx=parser.Results.dx;
    		p_left=norm_fun_cdf(x-dx,mu,v,fun,'lower',varargin{:});
    		p_right=norm_fun_cdf(x+dx,mu,v,fun,'lower',varargin{:});
    		f=max((p_right-p_left)/(2*dx),0);
        end
    elseif isstruct(fun)
        [w,k,lambda,s,m]=norm_quad_to_gx2_params(mu,v,fun);
        f=gx2pdf(x,w,k,lambda,s,m,varargin{:});
    end
end