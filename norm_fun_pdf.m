function f=norm_fun_pdf(x,mu,v,fun,varargin)

% parse inputs
parser = inputParser;
addRequired(parser,'x',@isnumeric);
addRequired(parser,'mu',@isnumeric);
addRequired(parser,'v',@isnumeric);
addRequired(parser,'fun',@(x) isa(x,'function_handle'));
addParameter(parser,'fun_span',3);
addParameter(parser,'fun_resol',100);
addParameter(parser,'AbsTol',1e-10);
addParameter(parser,'RelTol',1e-2);
addParameter(parser,'dx',1e-3,@(x) isreal(x) && isscalar(x) && (x>=0));

parse(parser,x,mu,v,fun,varargin{:});
dx=parser.Results.dx;
p_left=norm_fun_cdf(x-dx,mu,v,fun,varargin{:});
p_right=norm_fun_cdf(x+dx,mu,v,fun,varargin{:});
f=max((p_right-p_left)/(2*dx),0);

end