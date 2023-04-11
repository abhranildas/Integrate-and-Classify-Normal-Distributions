function [init_sign,x,samp_correct]=combine_ray_trace_doms(domlist,op,n,varargin)
% Combine multiple domains, defined as ray-trace functions, into one
% ray-trace domain, by intersection or union.
% Given an array of rays n from a specified origin, this
% function returns the initial signs (indicating if the beginning of the
% rays, i.e. -∞, are inside the combined domain) and the crossing points. The
% third output argument, samp_correct, treats n as an array of sample
% points, and returns whether the samples lie within the combined domain. This
% is used when classifying samples using this domain.
%
% Abhranil Das <abhranil.das@utexas.edu>
% Center for Perceptual Systems, University of Texas at Austin
% If you use this code, please cite:
% <a href="matlab:web('https://arxiv.org/abs/2012.14331')"
% >A method to integrate and classify normal distributions</a>.
%
% Inputs:
% domlist       cell array of ray-traced domain functions
% op            combination operation. 'and'=intersection, 'or'=union.
% n             array of ray directions or sample points in each column
% orig          column vector of the origin of rays
%
% Outputs:
% init_sign     initial signs of the combined domain along each ray
% x             cell array of crossing points of the combined domain along the
%               ray(s). This is a cell array because there may be
%               varying numbers of crossing points in each direction.
% samp_correct  this treats n as an array of samples, and returns if
%               they lie within the combined domain. If this is called, the
%               first two outputs are not returned.

% parse inputs
dim=size(n,1);
parser=inputParser;
addRequired(parser,'domlist');
addRequired(parser,'op');
addRequired(parser,'n',@isnumeric);
addParameter(parser,'orig',zeros(dim,1),@isnumeric);

parse(parser,domlist,op,n,varargin{:});
orig=parser.Results.orig;

n_doms=length(domlist);
n_dirs=size(n,2);

if nargout==3
    all_samp_correct=nan(n_doms,n_dirs);
    for i=1:n_doms
        [~,~,samp_correct_each]=domlist{i}(n,[]);
        all_samp_correct(i,:)=samp_correct_each;
    end
    if strcmpi(op,'or')
        samp_correct=any(all_samp_correct);
    elseif strcmpi(op,'and')
        samp_correct=all(all_samp_correct);
    end
    init_sign=[];
    x=[];

else
    % compile distances and signs to all boundaries
    all_init_sign=nan(n_doms,n_dirs); all_x=cell(1,n_dirs);
    for i=1:n_doms
        [init_sign_each,x_each]=domlist{i}(n,orig);
        all_x=arrayfun(@(x,y) [x{:};y], all_x,x_each,'un',0);
        all_init_sign(i,:)=init_sign_each;
    end

    [init_sign,x]=cellfun(@(all_init_sign_ray,all_x_ray) combine_ray_trace_doms_ray(all_init_sign_ray,all_x_ray,op),num2cell(all_init_sign,1),all_x,'un',0);
    init_sign=cell2mat(init_sign);
end