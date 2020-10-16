function p=norm_fun_cdf(x,mu,v,fun,varargin)

% parse inputs
parser=inputParser;
parser.KeepUnmatched=true;
addRequired(parser,'x',@isnumeric);
addRequired(parser,'mu',@isnumeric);
addRequired(parser,'v',@isnumeric);
addRequired(parser,'fun',@(x) isa(x,'function_handle'));
addOptional(parser,'side','lower',@(x) strcmpi(x,'lower') || strcmpi(x,'upper') );
addParameter(parser,'fun_span',3);
addParameter(parser,'fun_resol',100);
addParameter(parser,'AbsTol',1e-10);
addParameter(parser,'RelTol',1e-2);

parse(parser,x,mu,v,fun,varargin{:});
side=parser.Results.side;

[pc,p]=arrayfun(@(fun_level) integrate_normal(mu,v,fun,'reg_type','fun','fun_level',fun_level,'bPlot',false,varargin{:}), x);

if strcmpi(side,'upper')
    p=pc;
end
end
