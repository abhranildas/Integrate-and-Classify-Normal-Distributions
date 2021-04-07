close all;
clear all;
%rng(0);  % fix random seed

% First, generate an arbitrary R3 quadratic that looks kinda like a spherical obstacle
r = 2.5;
c = [0; 2; 1];
Q = -eye(3);
q = -Q*c;
q0 = 0.5*(r*r + c'*Q*c);

% Test reasonability
assert( 0.5*c'*Q*c + q'*c + q0 > 0 );   % should be positive at center
e = rand(3,1); e = e / norm(e);  % rand uvec
u = c + e*r;  % random point on surface
assert( abs( 0.5*u'*Q*u + q'*u + q0 ) < 1e-8 ); % should be zero on surface
u = c + 2*e*r;  % random point outside
assert( 0.5*u'*Q*u + q'*u + q0 < 0 ); % should be negative outside sphere

% Generate a random state distribution in position / velocity space
T = 0.2;  % length of time interval
mu_p = c + 1.4*r*e;  % position slightly outside sphere
mu_v = -0.5*r/T*e;  % velocity pointing into sphere

mu = [mu_p; mu_v];
P = 0.5*rand(6,6); P = P*P';  % rand PSD
mu_p
mu_v

% Integration domain is intersection of two constraints: q( p(0) ) < 0 AND q( p(T) ) > 0.
% Will specify this domain two ways: via lazy 'fun' and via explicit 'ray_scan' (more accurate).
% Concept is taken from following paper
%   Frey, Kristoffer M., Ted J. Steiner, and Jonathan P. How.
%   "Collision Probabilities for Continuous-Time Systems Without Sampling."
%   Robotics: Science and Systems (RSS), Corvallis, Oregon (2020).

% Wrap functions below:
dom = @(p1_,p2_,p3_,v1_,v2_,v3_) interval_collision_detect([p1_; p2_; p3_], [v1_; v2_; v3_], T, Q, q, q0);
dom_rays = @(nu_, mu_, P_) interval_collision_detect_ray_scan_normalized(nu_, mu_, P_, T, Q, q, q0);

s = interval_collision_detect(mu_p, mu_v, T, Q, q, q0);
fprintf('Domain value at mean: %0.3f -- positive if mean state is currently "safe" but headed for collision by time T.\n', s);
P
F0 = integrate_normal(mu, P, dom, 'dom_type', 'fun', 'mc_samples',1e4,'plotmode', 0);
F1 = integrate_normal(mu, P, dom_rays, 'dom_type', 'ray_scan', 'plotmode', 0);
fprintf('Probability of \"new\" collision: exact ray scan %0.3f vs numerically-sampled %0.3f.\n', F0, F1);


% Straightforward specification of the integration domain in R^6 (position + velocity).
function [g, g1, g2] = interval_collision_detect(p, v, T, Q, q, q0)
  assert(size(p,1) == 3);
  assert(size(v,1) == 3);
  assert(isequal(size(Q), [3,3]));
  assert(isequal(size(q), [3,1]));
  assert(isequal(size(q0), [1,1]));

  % Compute y0, y1
  p1 = p + T*v;
  y0 = 0.5*diag(p'*Q*p)' + q'*p + q0;
  y1 = 0.5*diag(p1'*Q*p1)' + q'*p1 + q0;

  g1 = -y0;
  g2 = y1;
  g = min(g1, g2);
end

% Wraps ray-scanning specification: takes rays sampled from normalized space, and warps to our "user space".
function [init_sign, Z] = interval_collision_detect_ray_scan_normalized(nu, mu, V, T, Q, q, q0)
  [n,N] = size(nu);
  assert(isequal( size(mu), [n,1] ));
  assert(isequal( size(V), [n,n] ));
  assert(isequal( size(T), [1,1] ));

  % Warp normals into domain space via sqrt Covariance.
  % Note: warped normals eta are *un*-normalized.
  L = chol(V);
  eta = L*nu;
  assert(isequal( size(eta), [n,N] ));

  % Do ray scan in domain space
  % No output rescaling needed for now...
  [init_sign, Z] = interval_collision_detect_ray_scan(eta, mu, T, Q, q, q0);
end

% "Ideal" user-level ray scan interface, all specificed in original space.
function [init_sign, Z] = interval_collision_detect_ray_scan(eta, x0, T, Q, q, q0)
  N = size(eta,2);
  assert(isequal(size(eta), [6,N]));
  peta = eta(1:3,:);
  assert(isequal(size(x0), [6,1]));
  p0 = x0(1:3);
  v0 = x0(4:6);
  assert(isequal(size(Q), [3,3]));
  assert(isequal(size(q), [3,1]));
  assert(isequal(size(q0), [1,1]));

  % Compute centered coefficients of g1, g2 with respect to ray length lambda
  g1 = [-0.5*diag(peta'*Q*peta)';
        -(q' + p0'*Q)*peta;
        repmat(-q0 - q'*p0 - 0.5*p0'*Q*p0, 1,N)];
  assert(isequal( size(g1), [3,N] ));

  Qbar = [Q, T*Q; T*Q, T*T*Q];
  qbar = [q; T*q];
  g2 = [ 0.5*diag(eta'*Qbar*eta)';
        (qbar' + x0'*Qbar)*eta;
        repmat(q0 + qbar'*x0 + 0.5*x0'*Qbar*x0, 1,N)];

  % TMP: Test compare with explicit interval_collision_detect()
  % Expect that p1(lambda) = g1(x0 + lambda*n)
  [~, g10, g20] = interval_collision_detect(p0, v0, T, Q, q, q0);
  assert( all( abs(g10 - g1(3,:)) < 1e-7 ) );
  assert( all( abs(g20 - g2(3,:)) < 1e-7 ) );
  for m = 1:N
    lambda = rand();
    p = p0 + lambda*eta(1:3,m);
    v = v0 + lambda*eta(4:6,m);

    [~, g1x, g2x] = interval_collision_detect(p, v, T, Q, q, q0);
    assert( abs(g1x - polyval(g1(:,m)', lambda)) < 1e-7 );
    assert( abs(g2x - polyval(g2(:,m)', lambda)) < 1e-7 );
  end

  [init_sign, Z] = ray_scan_joint_quad(g1, g2);
end

% Utility takes two *scalar*-valued quadratics and determines crossings of "joint" constraint min(g1, g2) > 0.
function [init_sign, Z] = ray_scan_joint_quad(g1, g2)
  N = size(g1, 2);
  assert( isequal(size(g1), [3,N]) );
  assert( isequal(size(g2), [3,N]) );

  % Unpack
  g1_a = g1(1,:);
  g1_b = g1(2,:);
  g1_c = g1(3,:);
  g2_a = g2(1,:);
  g2_b = g2(2,:);
  g2_c = g2(3,:);

  % Init sign is defined as sign at lambda = -inf! Not at origin!
  % Because g1, g2 quadratic, sign at -inf can be determined from coefficients,
  % starting with a, then -b, then c.
  s1 = sign(g1_a); s1(~s1) = -sign(g1_b(~s1)); s1(~s1) = sign(g1_c(~s1));
  s2 = sign(g2_a); s2(~s2) = -sign(g2_b(~s2)); s2(~s2) = sign(g2_c(~s2));
  init_sign = min(s1, s2);

  % Helper computes crossings (subset of roots) for single quadratic p
  function [z] = crossings(p)
    assert(isequal(size(p), [1,3]));
    z = roots(p);  % get roots
    z = z(~imag(z));  % filter out non-real roots

    % Compute slope at remaining roots and remove "grazing" points
    dp = p(1:2);
    dp(1) = 2*dp(1);
    slope = polyval(dp, z);
    z = z(slope ~= 0);
  end

  % Helper returns sign crossings for "joint" function
  function [z] = joint_crossings(g1_a, g2_a, g1_b, g2_b, g1_c, g2_c)
    % Collapse poly coefficients
    p1 = [g1_a, g1_b, g1_c];
    p2 = [g2_a, g2_b, g2_c];

    % Get crossing points and combine
    z1 = crossings(p1);
    z2 = crossings(p2);
    z = sort([z1; z2]);

    if ~isempty(z)
      % Generate test points in each interval
      L = z(1:end-1) + 0.5*diff(z);  % interior points
      x = [z(1)-1.0; L; z(end)+1.0];
      y = min( polyval(p1, x), polyval(p2, x) );
      z = z(diff(y > 0) ~= 0);  % indicates global crossings
    end

    % And return z as a row
    z = z';
  end

  % And compute crossings for each scalar poly pair
  Z = arrayfun(@joint_crossings, g1_a, g2_a, g1_b, g2_b, g1_c, g2_c, 'UniformOutput', false);
end
