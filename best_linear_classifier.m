function dom = best_linear_classifier(mu_1, v_1, mu_2, v_2, prior_1, vals)
% Ensure means are column vectors
mu_1 = mu_1(:); mu_2 = mu_2(:);

% Extract priors and calculate value differentials
p1 = prior_1;
p2 = 1 - prior_1;
dV1 = vals(1,1) - vals(1,2); % Value of predicting C1 over C2 when true is C1
dV2 = vals(2,2) - vals(2,1); % Value of predicting C2 over C1 when true is C2

% Constant K for shifting the threshold
K = (p2 * dV2) / (p1 * dV1);

% Find optimal t to minimize expected risk
t_opt = fminbnd(@calc_err, 0, 1, optimset('Display', 'off'));

% Get optimal projection (eta) and boundary (c) using the optimal t
[~, eta, c] = calc_err(t_opt);

% Orient the decision rule: a'*x + b > 0 for Class 1
dir = sign(eta'*mu_1 - eta'*mu_2);
if dir == 0; dir = 1; end % Fallback for exact overlap

a = dir * eta;
b = -dir * c;

dom.q2 = zeros(size(v_1));
dom.q1 = a;
dom.q0 = b;

% --- Nested helper function for optimization ---
    function [err, eta, c] = calc_err(t)
        % 1. Parametrize and normalize
        v = (t * v_1 + (1 - t) * v_2) \ (mu_2 - mu_1);
        eta = v / norm(v);

        % 2. Project distributions
        m1 = eta'*mu_1; v1 = eta'*v_1*eta;
        m2 = eta'*mu_2; v2 = eta'*v_2*eta;

        dir_t = sign(m1 - m2);
        if dir_t == 0; dir_t = 1; end

        % 3. Locate intersection boundary
        if abs(v1 - v2) < 1e-9
            B = 2*(m2*v1 - m1*v2);
            C = m1^2*v2 - m2^2*v1 + v1*v2*log(v1/v2) + 2*v1*v2*log(K);
            c = -C / B;

            % Direction-aware error (no abs to handle c outside means)
            err = p1*dV1*normcdf(dir_t*(c-m1)/sqrt(v1)) + p2*dV2*normcdf(dir_t*(m2-c)/sqrt(v2));
        else
            A = v2 - v1;
            B = 2*(m2*v1 - m1*v2);
            % Incorporate Bayesian threshold shift into C
            C = m1^2*v2 - m2^2*v1 + v1*v2*log(v1/v2) + 2*v1*v2*log(K);
            disc = B^2 - 4*A*C;

            if disc < 0
                % No intersection; one class's value heavily dominates the other everywhere
                err1 = p1*dV1; % Risk if always predicting Class 2
                err2 = p2*dV2; % Risk if always predicting Class 1
                if err1 < err2
                    c = -Inf * dir_t; err = err1;
                else
                    c = Inf * dir_t; err = err2;
                end
            else
                % Quadratic roots
                R = [(-B + sqrt(disc))/(2*A), (-B - sqrt(disc))/(2*A)];

                % Evaluate risk for both roots; strict between-means assumption is dropped
                err_R1 = p1*dV1*normcdf(dir_t*(R(1)-m1)/sqrt(v1)) + p2*dV2*normcdf(dir_t*(m2-R(1))/sqrt(v2));
                err_R2 = p1*dV1*normcdf(dir_t*(R(2)-m1)/sqrt(v1)) + p2*dV2*normcdf(dir_t*(m2-R(2))/sqrt(v2));

                if err_R1 < err_R2
                    c = R(1); err = err_R1;
                else
                    c = R(2); err = err_R2;
                end
            end
        end
    end
end