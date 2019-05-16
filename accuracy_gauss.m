function [mass,bd_pts,bd_sign]=accuracy_gauss(mu_a,vcov_a,mu_b,vcov_b,p_a,custom_bd_coeffs)

dim=length(mu_a);

% Cholesky decomposition of distribution a:
C_a=chol(vcov_a,'lower');

% if custom boundary is provided,
if ~isempty(custom_bd_coeffs)
    a2_c=custom_bd_coeffs.a2;
    a1_c=custom_bd_coeffs.a1;
    a0_c=custom_bd_coeffs.a0;
    
    % transform custom boundary
    a2=C_a'*a2_c*C_a;
    a1=C_a'*(2*a2_c*mu_a+a1_c);
    a0=mu_a'*a2_c*mu_a+a1_c'*mu_a+a0_c;
else
    % transformed distribution b:
    mu=C_a\(mu_b-mu_a);    
    S=C_a\vcov_b*inv(C_a)';
    
    %vx=sig(1,1); vy=sig(2,2); c=sig(1,2); D=det(sig);
    
    % boundary coefficients:
    a2=inv(S)-eye(dim);
    a1=-2*(S\mu);
    a0=mu'/S*mu+log((p_a/(1-p_a))^2*det(S));
end

% Create grid of unit vectors
if dim==1
    n_list=[1];
elseif dim==2
    dth=1e-3;
    th=-pi/2:dth:pi/2;
    n_list=[cos(th);sin(th)];    
elseif dim==3    
    dth=5e-2; dph=5e-2;
    th=0:dth:pi/2;
    ph=0:dph:2*pi;
    n_list=[];
    for i=1:length(th)
        for j=1:length(ph)
            n_list=[n_list,[sin(th(i))*cos(ph(j)); sin(th(i))*sin(ph(j)); cos(th(i))]];
        end
    end
%     n_list=[sin(th).*cos(ph);sin(th).*sin(ph);cos(th)];
end

% Integrate
bd_pts_white=[];
bd_sign=[];

mass=0;
if a0>0
    mass=1;
end

m=0;
for i=1:length(n_list)
    n=n_list(:,i);
    q2=n'*a2*n;
    q1=a1'*n;
    
    %     q=  [vy+vx*k^2+2*c*k-D*(1+k^2),...
    %         2*c*(k*mx+my)-mx*vy-k*my*vx,...
    %         vx*my^2+vy*mx^2+2*c*mx*my+log(D^D)];
    
    rs=roots([q2 q1 a0])'; % quadratic coefficients
    rs=rs(~imag(rs)); % only real roots
    for r=rs
        coords=n*r;
        bd_pts_white=[bd_pts_white,coords];
        r_sign=sign(2*q2*r^2+q1*r);
        bd_sign=[bd_sign, r_sign];
        if dim==1            
            dm=r_sign*normcdf(abs(r),'upper');
        elseif dim==2
            dm=r_sign*exp(-r^2/2);
        elseif dim==3
            sinth=sqrt(1-n(3)^2);
            dm=r_sign * (sqrt(pi)/2 * (1-erf(abs(r)/sqrt(2))) + abs(r)*exp(-abs(r)^2/2)/sqrt(2))*sinth;
        end
        m=m+sum(dm);
    end
end

% factors
if dim==2
    m=m*dth/(2*pi);
elseif dim==3
    m=m*dth*dph/(2*pi^(3/2));
end

mass=mass+m;

% transform boundary points to raw space:
if numel(bd_pts_white)
    bd_pts=C_a*bd_pts_white+mu_a;
else
    bd_pts=[];
end