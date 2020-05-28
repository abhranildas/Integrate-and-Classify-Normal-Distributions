function f=gx2pdf(x,lambda,m,delta,c,varargin)
% Returns the PDF of a generalized chi-squared (a weighted sum of
% non-central chi-squares).

% Syntax:
% f=gx2pdf(x,lambda,m,delta,c)
% f=gx2pdf(x,lambda,m,delta,c,'dx',1e-1)
% f=gx2pdf(x,lambda,m,delta,c,'AbsTol',0,'RelTol',1e-7)
% f=gx2pdf(x,lambda,m,delta,c,'approx','tail')

% Example:
% f=gx2pdf(25,[1 -5 2],[1 2 3],[2 3 7],0)

% Inputs:
% x         point at which to evaluate the PDF
% lambda    row vector of coefficients of the non-central chi-squares
% m         row vector of degrees of freedom of the non-central chi-squares
% delta     row vector of non-centrality paramaters (sum of squares of
%           means) of the non-central chi-squares
% c         constant term
% dx        fineness for numerically differentiating the CDF to compute PDF
% 'AbsTol'  absolute error tolerance for the CDF computation
% 'RelTol'  relative error tolerance for the CDF computation
%           The absolute OR the relative tolerance is satisfied.
% 'approx'  set to 'tail' for Pearson's approximation of the CDF tail.
%           Works best for the upper (lower) tail when all
%           lambda are positive (negative).

% Output:
% f         computed PDF

% Author:
% Abhranil Das <abhranil.das@utexas.edu>
% Center for Perceptual Systems, University of Texas at Austin

% If you use this code, you may cite:
% A new method to compute classification error
% jov.arvojournals.org/article.aspx?articleid=2750251

parser = inputParser;
addRequired(parser,'x',@(x) isreal(x) && isscalar(x));
addRequired(parser,'lambda',@(x) isreal(x) && isrow(x));
addRequired(parser,'m',@(x) isreal(x) && isrow(x));
addRequired(parser,'delta',@(x) isreal(x) && isrow(x));
addRequired(parser,'c',@(x) isreal(x) && isscalar(x));
addParameter(parser,'dx',1e-10,@(x) isreal(x) && isscalar(x) && (x>=0));
addParameter(parser,'AbsTol',1e-10,@(x) isreal(x) && isscalar(x) && (x>=0));
addParameter(parser,'RelTol',1e-6,@(x) isreal(x) && isscalar(x) && (x>=0));
addParameter(parser,'approx','none',@(x) strcmpi(x,'none') || strcmpi(x,'tail'));

parse(parser,x,lambda,m,delta,c,varargin{:});

dx=parser.Results.dx;
if any(strcmp(varargin,'dx'))
    removeIndex=strcmp(varargin(:,1),'dx');
    varargin(removeIndex,:)=[];
end

p_right=gx2cdf_imhof(x+dx,lambda,m,delta,c,varargin{:});
p_left=gx2cdf_imhof(x-dx,lambda,m,delta,c,varargin{:});

f=(p_right-p_left)/(2*dx);

end