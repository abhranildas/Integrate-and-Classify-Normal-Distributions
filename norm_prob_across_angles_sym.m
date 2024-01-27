function [p_angles,bd_pts_angles]=norm_prob_across_angles_sym(mu,v,dom,theta,varargin)
    theta=double(theta);
    [p_angles,bd_pts_angles]=norm_prob_across_angles(mu,v,dom,'theta',theta,varargin{:});