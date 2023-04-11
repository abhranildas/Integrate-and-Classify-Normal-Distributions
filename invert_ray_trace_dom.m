function [init_sign,x,samp_correct]=invert_ray_trace_dom(dom,n,varargin)
% Invert a domain defined as ray-trace functions.
% Given an array of rays n from a specified origin, this
% function returns the initial signs (indicating if the beginning of the
% rays, i.e. -∞, are inside the domain) and the crossing points. The
% third output argument, samp_correct, treats n as an array of sample
% points, and returns whether the samples lie within the domain. This
% is used when classifying samples using this domain.
%
% Abhranil Das <abhranil.das@utexas.edu>
% Center for Perceptual Systems, University of Texas at Austin
% If you use this code, please cite:
% <a href="matlab:web('https://arxiv.org/abs/2012.14331')"
% >A method to integrate and classify normal distributions</a>.
%
% Inputs:
% dom           ray-traced domain functions
% n             array of ray directions or sample points in each column
% orig          column vector of the origin of rays
%
% Outputs:
% init_sign     initial signs of the domain along each ray
% x             cell array of crossing points of the domain along the
%               ray(s). This is a cell array because there may be
%               varying numbers of crossing points in each direction.
% samp_correct  this treats n as an array of samples, and returns if
%               they lie within the domain. If this is called, the
%               first two outputs are not returned.

% parse inputs
dim=size(n,1);
parser = inputParser;
addRequired(parser,'dom');
addRequired(parser,'n',@isnumeric);
addParameter(parser,'orig',zeros(dim,1),@isnumeric);

parse(parser,dom,n,varargin{:});
orig=parser.Results.orig;

if nargout==3
    [~,~,samp_correct]=dom(n,[]);
    samp_correct=~samp_correct;
    init_sign=[];
    x=[];    
else    
    [init_sign,x]=dom(n,orig);
    init_sign=-init_sign;
end