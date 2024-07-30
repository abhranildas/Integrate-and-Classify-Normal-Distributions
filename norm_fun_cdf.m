function p=norm_fun_cdf(x,mu,v,fun,varargin)
% NORM_FUN_CDF Compute the cdf of a scalar function of a (multi)normal
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
% norm_fun_cdf(0.5,mu,v,fun)
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
% Optional positional input:
% 'upper'       more accurate complementary cdf
%
% Optional name-value inputs:
% fun_span      trace radius (in Mahalanobis distance) for function.
%               Default=5.
% fun_resol     resolution of tracing (finding roots) of function.
%               Default=100.
% AbsTol        absolute tolerance for the output
% RelTol        relative tolerance for the output
%               The absolute OR the relative tolerance will be satisfied.
%
% Outputs:
% p             cdf
%
% See also:
% <a href="matlab:open(strcat(fileparts(which('integrate_normal')),filesep,'doc',filesep,'GettingStarted.mlx'))">Interactive demos</a>
% norm_fun_pdf
% norm_fun_inv

% parse inputs
parser=inputParser;
parser.KeepUnmatched=true;
addRequired(parser,'x',@isnumeric);
addRequired(parser,'mu',@isnumeric);
addRequired(parser,'v',@isnumeric);
addRequired(parser,'fun',@(x) isstruct(x)|| isa(x,'function_handle'));
addOptional(parser,'side','lower',@(x) strcmpi(x,'lower') || strcmpi(x,'upper') );
addParameter(parser,'AbsTol',1e-10);
addParameter(parser,'RelTol',1e-2);

parse(parser,x,mu,v,fun,varargin{:});

if isa(fun,'function_handle')
    % if strcmpi(side,'lower')
    if ismember('side',parser.UsingDefaults)
        p=arrayfun(@(fun_level) integrate_normal(mu,v,fun,'lower',varargin{:},'dom_type','fun','fun_level',fun_level,'plotmode',false), x);
    else
        p=arrayfun(@(fun_level) integrate_normal(mu,v,fun,varargin{:},'dom_type','fun','fun_level',fun_level,'plotmode',false), x);
    end        
    % elseif strcmpi(side,'upper')
        % p=arrayfun(@(fun_level) integrate_normal(mu,v,fun,varargin{:},'dom_type','fun','fun_level',fun_level,'dom_side','above','plotmode',false), x);
    % end
elseif isstruct(fun)
    [w,k,lambda,s,m]=norm_quad_to_gx2_params(mu,v,fun);
    p=gx2cdf(x,w,k,lambda,s,m,varargin{:});
end
end
