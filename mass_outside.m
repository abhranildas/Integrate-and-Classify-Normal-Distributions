function [m,dec_bd_raw,dec_bd_sign]=mass_outside(mu_a,sig_a,mu_b,sig_b,p_a)

d=length(mu_a);

if ~exist('p_a','var')
    p_a=.5;
end

% Cholesky decomposition of distribution a:
C_a=chol(sig_a,'lower');

% transformed distribution b:
mu=C_a\(mu_b-mu_a);
%mu=(mu_b-mu_a);
%mx=mu(1); my=mu(2);

S=C_a\sig_b*inv(C_a)';

%vx=sig(1,1); vy=sig(2,2); c=sig(1,2); D=det(sig);

% matrix coefficients of decision-boundary quadratic:
a2=eye(d)-inv(S);
a1=2*mu'/S;
a0=-mu'/S*mu-log((p_a/(1-p_a))^2*det(S));

% Create grid of unit vectors
if d==2
    dth=1e-3;
    th=-pi/2:dth:pi/2;
    n_list=[cos(th);sin(th)];    
elseif d==3    
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
dec_bd=[];
dec_bd_sign=[];
m=0;
for i=1:length(n_list)
    n=n_list(:,i);
    q2=n'*a2*n;
    q1=a1*n;
    
    %     q=  [vy+vx*k^2+2*c*k-D*(1+k^2),...
    %         2*c*(k*mx+my)-mx*vy-k*my*vx,...
    %         vx*my^2+vy*mx^2+2*c*mx*my+log(D^D)];
    
    r=roots([q2 q1 a0]); % quadratic coefficients
    if ~isempty(r)&&isreal(r) % if there is a solution
        coords=n*r';
        dec_bd=[dec_bd,coords];
        signs=sign(2*q2*r.^2+q1*r); % 0 to pt, or pt to infty
        dec_bd_sign=[dec_bd_sign, signs'];
        if d==2            
            m=m+sum(signs.*(exp(-r.^2/2)-.5)+.5); 
        elseif d==3
            sinth=sqrt(1-n(3)^2);
            r=abs(r);
            dm=sum(sqrt(pi)/4+...
                signs.*(r.*exp(-r.^2/2)/sqrt(2)-sqrt(pi)/2*erf(r/sqrt(2))+sqrt(pi)/4))*sinth;
            m=m+dm;
        end
    end
end

if d==2
    m=m*dth/(2*pi);
elseif d==3
    m=m*dth*dph/(2*pi^(3/2));
end

% reverse-transform decision bound:
dec_bd_raw=C_a*dec_bd+repmat(mu_a,[1 size(dec_bd,2)]);