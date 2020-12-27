function p=norm_vfun_cdf(f,mu,v,vfun,varargin)

% f is a matrix, where each row is a point to evaluate the cdf, like
% mvncdf
% parse inputs
parser=inputParser;
parser.KeepUnmatched=true;
addRequired(parser,'f',@isnumeric);
addRequired(parser,'mu',@isnumeric);
addRequired(parser,'v',@isnumeric);
addRequired(parser,'vfun',@(x) isa(x,'function_handle'));
addOptional(parser,'side','lower',@(x) strcmpi(x,'lower') || strcmpi(x,'upper') );
addParameter(parser,'fun_span',3);
addParameter(parser,'fun_resol',100);
addParameter(parser,'AbsTol',1e-10);
addParameter(parser,'RelTol',1e-2);

parse(parser,f,mu,v,vfun,varargin{:});
% side=parser.Results.side;

dvfun=@(x,y,f_pt) arrayfun(@(x_pt,y_pt) max(vfun(x_pt,y_pt)-f_pt),x,y); % scalar decision variable function from vector function


p=arrayfun(@(row_idx) norm_fun_cdf(0,mu,v,@(x,y) dvfun(x,y,f(row_idx,:)),varargin{:}),(1:size(f,1))');

% [pc,p]=integrate_normal(mu,v,dvfun,'dom_type','fun','plotmode',false,varargin{:});

% arrayfun(@(x,y) norm_vfun_cdf([x y],mu,v,f3),X,Y);

end
