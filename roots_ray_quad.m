function x_ray=roots_ray_quad(q2_pt,q1_pt,q0_pt)
    % root(s) along a ray through quad domain
    x_ray=sort(roots([q2_pt q1_pt q0_pt])).';
    x_ray=x_ray(~imag(x_ray)); % only real roots

    % remove any roots that are tangents. Only crossing points
    slope=2*q2_pt*x_ray+q1_pt;
    x_ray=x_ray(slope~=0);
end