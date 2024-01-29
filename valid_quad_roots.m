function r=valid_quad_roots(q2,q1,q0)
    r=sort(roots([q2 q1 q0])).';
    r=r(~imag(r)); % only real roots

    % remove any roots that are tangents. Only actual crossing points
    slope=2*q2*r+q1;
    r=r(slope~=0);
end