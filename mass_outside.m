function [m,dec_bd_raw,dec_bd_sign]=mass_outside(mu_a,sig_a,mu_b,sig_b,dth)

% Cholesky decomposition of first distribution:
C_a=chol(sig_a,'lower');

% Transformed second distribution:
mu=inv(C_a)*(mu_b-mu_a);
%mu=(mu_b-mu_a);
%mx=mu(1); my=mu(2);

sig=inv(C_a)*sig_b*inv(C_a)';

%vx=sig(1,1); vy=sig(2,2); c=sig(1,2); D=det(sig);

% Terms for decision-bound quadratic:
a2_mat=inv(sig)-eye(2);
a1_mat=-mu'*(inv(sig)+(inv(sig))');
a0=mu'/sig*mu+log(det(sig));

if ~exist('dth','var')
    dth=1e-2;
end

m=0;
th_grid=-pi/2:dth:pi/2;
dec_bd=[];%nan(2,length(th_grid));
dec_bd_sign=[];

for i=1:length(th_grid)
    th=th_grid(i);
    k=tan(th);   
    k_v=[1;tan(th)];
    a2=k_v'*a2_mat*k_v;
    a1=a1_mat*k_v;
    
    %     q=  [vy+vx*k^2+2*c*k-D*(1+k^2),...
    %         2*c*(k*mx+my)-mx*vy-k*my*vx,...
    %         vx*my^2+vy*mx^2+2*c*mx*my+log(D^D)];
    q=[a2 a1 a0];
    x=roots(q);
    if isreal(x) % if there is a solution
        y=k*x;
        dec_bd=[dec_bd,[x,y]'];
        signs=sign(2*a2*x.^2+a1*x);
        dec_bd_sign=[dec_bd_sign, signs'];
        r2=x.^2+y.^2;
        m=m+sum(signs.*(exp(-r2/2)-.5)+.5); % 0 to pt, or pt to infty
    end
end

m=m*dth/(2*pi);

% figure; hold on
% plot(z(1,:),z(2,:),'-b')
% plot(dist_b(1,:),dist_b(2,:),'-r')
% plot(dec_bd(1,:),dec_bd(2,:),'.g')
% axis image

% reverse-transform decision bound:
dec_bd_raw=C_a*dec_bd+repmat(mu_a,[1 size(dec_bd,2)]);