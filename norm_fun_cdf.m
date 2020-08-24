function p=norm_fun_cdf(f,mu,v,fun,varargin)
% Integrate a normal distribution over a specified region.
%
% How to use this command:
% See github readme at https://github.com/abhranildas/classify
%
% Credits:
%   Abhranil Das <abhranil.das@utexas.edu>
%	R Calen Walshe
%	Wilson S Geisler
%	Center for Perceptual Systems, University of Texas at Austin
% If you use this code, please cite:
%   A new method to compute classification error
%   https://jov.arvojournals.org/article.aspx?articleid=2750251

% parse inputs
parser = inputParser;
addRequired(parser,'f',@isnumeric);
addRequired(parser,'mu',@isnumeric);
addRequired(parser,'v',@isnumeric);
addRequired(parser,'fun',@(x) isstruct(x)|| isa(x,'function_handle'));
addOptional(parser,'side','lower',@(x) strcmpi(x,'lower') || strcmpi(x,'upper') );
addParameter(parser,'fun_span',3);
addParameter(parser,'fun_resol',100);
addParameter(parser,'AbsTol',1e-10);
addParameter(parser,'RelTol',1e-2);
addParameter(parser,'bPlot',0);

parse(parser,f,mu,v,fun,varargin{:});
fun=parser.Results.fun;
fun_span=parser.Results.fun_span;
fun_resol=parser.Results.fun_resol;
AbsTol=parser.Results.AbsTol;
RelTol=parser.Results.RelTol;
bPlot=parser.Results.bPlot;

if isa(fun,'function_handle')
    reg_type='fun';
    reg=@(varargin) cdf_fun(fun,f,varargin);
end
[p,pc]=integrate_normal(mu,v,reg,'reg_type',reg_type)
