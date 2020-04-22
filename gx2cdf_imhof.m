function p=gx2cdf_imhof(x,lambda,m,delta,varargin)
% Returns the CDF of a generalized chi-squared (a weighted sum of
% non-central chi-squares), using Imhof's [1961] algorithm.

% Syntax:
% p=gx2cdf_imhof(x,lambda,m,delta)
% p=gx2cdf_imhof(x,lambda,m,delta,'upper')
% p=gx2cdf_imhof(x,lambda,m,delta,'AbsTol',0,'RelTol',1e-7)
% p=gx2cdf_imhof(x,lambda,m,delta,'upper','approx','tail')

% Example:
% p=gx2cdf_imhof(25,[1 -5 2],[1 2 3],[2 3 7])

% Inputs:
% x         point at which to evaluate the CDF
% lambda    row vector of coefficients of the non-central chi-squares
% m         row vector of degrees of freedom of the non-central chi-squares
% delta     row vector of non-centrality paramaters (sum of squares of
%           means) of the non-central chi-squares
% 'upper'   more accurate estimate of the complementary CDF when it's small
% 'AbsTol'  absolute error tolerance for the output.
% 'RelTol'  relative error tolerance for the output.
%           The absolute OR the relative tolerance is satisfied.
% 'approx'  set to 'tail' for Pearson's approximation of the tail
%           probabilities. Works best for the upper (lower) tail when all
%           lambda are positive (negative).

% Output:
% p         computed CDF

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
addOptional(parser,'side','lower',@(x) strcmpi(x,'lower') || strcmpi(x,'upper') );
addParameter(parser,'AbsTol',1e-10,@(x) isreal(x) && isscalar(x) && (x>=0));
addParameter(parser,'RelTol',1e-6,@(x) isreal(x) && isscalar(x) && (x>=0));
addParameter(parser,'approx','none',@(x) strcmpi(x,'none') || strcmpi(x,'tail'));

parse(parser,x,lambda,m,delta,varargin{:});
side=parser.Results.side;
approx=parser.Results.approx;

u=[]; % pre-allocate in static workspace

% define the integrand (lambda, m, delta must be column vectors here)
    function f=imhof_integrand(u,x,lambda,m,delta)
        theta=sum(m.*atan(lambda*u)+(delta.*(lambda*u))./(1+lambda.^2*u.^2),1)/2-u*x/2;
        rho=prod(((1+lambda.^2*u.^2).^(m/4)).*exp(((lambda.^2*u.^2).*delta)./(2*(1+lambda.^2*u.^2))),1);
        f=sin(theta)./(u.*rho);
    end

if strcmpi(approx,'tail') % compute tail approximations
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
    if any(strcmp(parser.UsingDefaults,'AbsTol')) && any(strcmp(parser.UsingDefaults,'RelTol'))
        imhof_integral=integral(@(u) imhof_integrand(u,x,lambda',m',delta'),0,inf);
        if strcmpi(side,'lower')
            p=0.5-imhof_integral/pi;
        elseif strcmpi(side,'upper')
            p=0.5+imhof_integral/pi;
        end
    else
        syms u
        imhof_integral=vpaintegral(@(u) imhof_integrand(u,x,lambda',m',delta'),u,0,inf,'AbsTol',parser.Results.AbsTol,'RelTol',parser.Results.RelTol,'MaxFunctionCalls',inf);
        
        if strcmpi(side,'lower')
            p=double(0.5-imhof_integral/pi);
        elseif strcmpi(side,'upper')
            p=double(0.5+imhof_integral/pi);
        end
    end
end

if p<0 || p>1
    warning("Inaccurate estimate of a small probability. Setting to NaN.")
    p=nan;
end

end