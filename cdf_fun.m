function out=cdf_fun(f,f0,x)
% parser = inputParser;
% parse(parser,f,x,varargin{:});
out=f(x{:})+f0;
