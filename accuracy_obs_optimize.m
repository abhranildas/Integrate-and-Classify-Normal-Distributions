function acc_obs=accuracy_obs_optimize(dim,x,obs_a,obs_b)

if dim==1
    a2=x(1);
    a1=x(2);
    a0=x(3);
elseif dim==2
    a2=[x(1) x(2); x(3) x(4)];
    a1=[x(5); x(6)];
    a0=x(7);
elseif dim==3
    a2=[x(1) x(2) x(3); x(4) x(5) x(6); x(7) x(8) x(9)];
    a1=[x(10); x(11); x(12)];
    a0=x(13);
end

acc_obs=-sum([dot(obs_a,obs_a*a2',2) + obs_a*a1 + a0 > 0;...
    dot(obs_b,obs_b*a2',2) + obs_b*a1 + a0 < 0]);