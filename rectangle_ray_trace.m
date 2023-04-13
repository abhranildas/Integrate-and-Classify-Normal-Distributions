function [init_sign,x] = rectangle_ray_trace(rect,n,orig)
% Return (hyper-)rectangle domains in ray-trace format.
% Given a ray from a specified origin in a vector direction n, this
% function returns the initial sign (indicating if the beginning of the
% ray, i.e. -∞, is inside the domain) and the crossing points.
%
% Shizhuang Wang <sz.wang@sjtu.edu.cn>
% School of Aeronautics and Astronautics, Shanghai Jiao Tong University
% using the method of 
% https://tavianator.com/2022/ray_box_boundary.html
%
% Abhranil Das <abhranil.das@utexas.edu>
% Center for Perceptual Systems, University of Texas at Austin
%
% If you use this code, please cite:
% <a href="matlab:web('https://arxiv.org/abs/2012.14331')"
% >A method to integrate and classify normal distributions</a>.
%
% Example:
% L = [-1 -1]; U = [1 1]; ray = [0 1;1 1]; orig=[0;0];
% [init_sign,cross_pts] = rectangle_ray_trace([L;U],ray,orig);
%
% Inputs:
% rect          matrix where the first and second row are lower and upper
%               limits of the (hyper-)rectangle
% n             array of ray directions or sample points in each column
% orig          column vector of the origin of the rays
%
% Outputs:
% init_sign     initial signs of the domain along each ray
% x             cell array of crossing points of the domain along the
%               ray(s). This is a cell array because there may be
%               varying numbers of crossing points in each direction.

%% Check Inputs
L=rect(1,:)'; U=rect(2,:)';
n_var=length(L);

% Check the sizes of {L,R,n,orig} don't match
if any([size(n,1),length(orig)] ~= n_var)
    error('Input dimensions do not match.');
end

%% Determine the initial sign
% (indicating if the beginning of the ray, i.e. -∞, is inside the domain)

n_ray = size(n,2);
init_sign = zeros(1,n_ray);

for i = 1 : n_ray
    init_pt = n(:,i) * (-Inf);
    if any(init_pt < L) || any(init_pt > U) % not in the domain
        init_sign(1,i) = -1;
    else                                    % in the domain
        init_sign(1,i) = 1;
    end
end

%% Determine the crossing points

% First, determine whether the origin is in the rectangle domain
flag = 1; % 1 : in the domain; 0: outside the domain
if any(orig<L) || any(orig>U) % outside the region
    flag = 0;
end

% Then, determine the crossing points
x = cell(1,n_ray);

% If the origin is inside the domain
if flag == 1

    for i = 1 : n_ray
        ray_i = n(:,i) / norm(n(:,i)); % normalization
        % calculate intersection with each axis-aligned hyperplane
        rou = zeros(n_var*2,1);
        for j = 1 : n_var
            if ray_i(j) ~= 0
                rou(j*2-1,1) =  (L(j) - orig(j)) / ray_i(j);
                rou(j*2,  1) =  (U(j) - orig(j)) / ray_i(j);
            else
                rou(j*2-1,1) = Inf;
                rou(j*2,  1) = Inf;
            end
        end
        % Two crossing points are in opposite directions
        rou_positive     = rou(rou>=0);
        rou_negative     = rou(rou< 0);
        cross_point(1,1) = max(rou_negative); % neareast
        cross_point(1,2) = min(rou_positive); % neareast
        x{1,i}   = cross_point;
    end

    % If the origin is outside the domain (efficient version)
    % reference: https://tavianator.com/2022/ray_box_boundary.html
else

    for i = 1 : n_ray

        ray_i = n(:,i) / norm(n(:,i)); % normalization
        ray_i_inv = 1./ray_i;
        rou = zeros(n_var,2);
        tmin_p = 0; tmax_p = Inf; % positive intersection
        tmin_n = 0; tmax_n = Inf; % negative intersection
        intersection = 0;
        for j = 1 : n_var
            rou(j,1) = (L(j) - orig(j)) * ray_i_inv(j); t1 = rou(j,1);
            rou(j,2) = (U(j) - orig(j)) * ray_i_inv(j); t2 = rou(j,2);
            tmin_pj = min(t1,t2); tmin_nj = min(-t1,-t2);
            tmax_pj = max(t1,t2); tmax_nj = max(-t1,-t2);
            tmin_p = max(tmin_p,tmin_pj);
            tmax_p = min(tmax_p,tmax_pj);
            tmin_n = max(tmin_n,tmin_nj);
            tmax_n = min(tmax_n,tmax_nj);
        end

        if tmin_p < tmax_p, intersection = 1;  end
        if tmin_n < tmax_n, intersection = -1; end

        if intersection == 1      % positive intersection
            x{1,i} = [tmin_p,tmax_p];
        elseif intersection == -1 % negative intersection
            x{1,i} = [-tmax_n,-tmin_n];
        else                      % no intersection
            x{1,i} = zeros(1,0);
        end
    end

end


