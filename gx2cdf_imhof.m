function p=gx2cdf_imhof(x,lambda,m,delta,varargin)
% Syntax:
% p=gx2cdf_imhof(x,lambda,m,delta)
% p=gx2cdf_imhof(x,lambda,m,delta,'tail')

% Description:
% Returns the CDF of a generalized chi-squared (a weighted sum of
% non-central chi-squares), using Imhof's [1961] algorithm.

% Example:
% [p,pc]=gx2cdf_imhof(25,[1 -5 2],[1 2 3],[2 3 7])

% Inputs:
% x         point at which to evaluate the CDF
% lambda    row vector of coefficients of the non-central chi-squares
% m         row vector of degrees of freedom of the non-central chi-squares
% delta     row vector of non-centrality paramaters (sum of squares of
%           means) of the non-central chi-squares

% Outputs:
% p         computed CDF

% Author:
% Abhranil Das <abhranil.das@utexas.edu>
% Center for Perceptual Systems, University of Texas at Austin

% If you use this code, you may cite:
% A new method to compute classification error
% jov.arvojournals.org/article.aspx?articleid=2750251

% define the integrand (lambda, m, delta must be column vectors here)

parser = inputParser;
addRequired(parser,'x',@(x) isnumeric(x) && isscalar(x));
addRequired(parser,'lambda',@isrow);
addRequired(parser,'m',@isrow);
addRequired(parser,'delta',@isrow);
addOptional(parser,'side','lower',@(x) strcmpi(x,'lower') || strcmpi(x,'upper') );
addParameter(parser,'estimate',[]);

parse(parser,x,lambda,m,delta,varargin{:});
side=parser.Results.side;
estimate=parser.Results.estimate;

if strcmpi(estimate,'tail') % compute tail approximations
    j=(1:3)';
    c=sum((lambda.^j).*(j.*delta+m),2);
    h=c(2)^3/c(3)^2;
    if c(3)>0
        y=(x-c(1))*sqrt(h/c(2))+h;
        if strcmpi(side,'lower')
            p=chi2cdf(y,h);
        elseif strcmpi(side,'upper')
            p=chi2cdf(y,h,'upper');
        end
    else
        c=sum(((-lambda).^j).*(j.*delta+m),2);
        y=(-x-c(1))*sqrt(h/c(2))+h;
        if strcmpi(side,'lower')
            p=chi2cdf(y,h,'upper');
        elseif strcmpi(side,'upper')
            p=chi2cdf(y,h);
        end
    end
    
else
    % compute the integral
    if isempty(estimate)
        imhof_integral=integral(@(u) imhof_integrand(u,x,lambda',m',delta'),0,inf);
        if strcmpi(side,'lower')
            p=0.5-imhof_integral/pi;
        elseif strcmpi(side,'upper')
            p=0.5+imhof_integral/pi;
        end
    else
        syms u
        imhof_integral=vpaintegral(@(u) imhof_integrand(u,x,lambda',m',delta'),u,0,inf,'RelTol',estimate,'AbsTol',0,'MaxFunctionCalls',inf);
    end
    
    if strcmpi(side,'lower')
        p=double(0.5-imhof_integral/pi);
    elseif strcmpi(side,'upper')
        p=double(0.5+imhof_integral/pi);
    end
    
    if p<1e-3 && isempty(estimate)
        warning("Tail probability is sometimes inaccurate, and might be improved by setting 'estimate' to a small tolerance (slow, accurate), or to 'tail' (fast, approximate, works best for upper tail with lambdas the same sign).");
    end
    
end

if p<0 || p>1
    warning("Inaccurate estimate of a small probability. Setting to NaN.")
    p=nan;
end

end