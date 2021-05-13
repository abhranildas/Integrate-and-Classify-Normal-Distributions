function [init_sign,x,samp_correct]=line_ray_trace(n,orig)
    % An example domain in 2d, y>k with k=1, in ray-trace format.
    % Given a ray from a specified origin in a vector direction n, this 
    % function returns the initial sign (indicating if the beginning of the
    % ray, i.e. -∞, is inside the domain) and the crossing points. The
    % third output argument, samp_correct, treats n as an array of sample
    % points, and returns whether the samples lie within the domain. This
    % is used when classifying samples using this domain.
    % For more details, see accompanying paper linked below.
    %
    % Abhranil Das <abhranil.das@utexas.edu>
    % Center for Perceptual Systems, University of Texas at Austin
    % If you use this code, please cite:
    % <a href="matlab:web('https://arxiv.org/abs/2012.14331')"
    % >A method to integrate and classify normal distributions</a>.
    %
    % Example:
    % n=[0 1 -1; 2 1 -1]; orig=[0;0];
    % [init_sign,x]=line_ray_trace(n,orig)
    % [~,~,samp_correct]=line_ray_trace(n,orig)
    %
    % Inputs:
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
    
    k=1;
    if nargout==3 % if the third output is being asked
        samp_correct=n(2,:)>k; % return if y coordinate of the samples is > k
        init_sign=[];
        x=[];
    else
        init_sign=sign(-n(2,:));
        n=n./vecnorm(n);
        x=num2cell((k-orig(2))./n(2,:));
    end