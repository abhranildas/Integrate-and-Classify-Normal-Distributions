function samp_bd_current=optimize_samp_bd(samp_1,samp_2,norm_bd,varargin)
parser=inputParser;
parser.KeepUnmatched=true;
addRequired(parser,'samp_1',@isnumeric);
addRequired(parser,'samp_2',@isnumeric);
addRequired(parser,'norm_bd',@isstruct);
addParameter(parser,'vals',eye(2), @(x) isnumeric(x) && ismatrix(x));
addParameter(parser,'samp_opt','svm');
addParameter(parser,'samp_opt_plot','step', @(x) strcmpi(x,'step') || strcmpi(x,'smooth') || x==false);
addParameter(parser,'plots',false, @islogical);

parse(parser,samp_1,samp_2,norm_bd,varargin{:});
samp_opt=parser.Results.samp_opt;
vals=parser.Results.vals;
plots=parser.Results.plots;

dim=size(samp_1,2);

% start from the normal boundary
norm_bd_flat=[norm_bd.q2(triu(true(size(norm_bd.q2)))); norm_bd.q1(:); norm_bd.q0];
samp_bd_current=norm_bd_flat;
samp_val_current=samp_value_flat(samp_bd_current,samp_1,samp_2,varargin{:});

if samp_opt
    if isscalar(samp_opt)
        obj_fun=@(x) -samp_value_flat(x,samp_1,samp_2,varargin{:}); % step objective function

        % if samp_opt==1, first try fminsearch (local optimization)
        [samp_bd_current,samp_val_best]=fminsearch(obj_fun,samp_bd_current,optimset('Display','notify','TolX',0,'TolFun',1e-2/(size(samp_1,1)+size(samp_2,1))));

        if samp_opt>1
            % then do global optimization, starting from the local fminsearch
            % step
            problem = createOptimProblem('fminunc', ...
                'x0',        samp_bd_current, ...
                'objective', obj_fun, ...
                'options',   optimset('Display','notify', ...
                'TolX', 0, ...
                'TolFun', 1e-2/(size(samp_1,1)+size(samp_2,1))));

            ms = MultiStart('UseParallel', true, 'StartPointsToRun', 'all','FunctionTolerance',0,'Display','final');

            [samp_bd_current, samp_val_best] = run(ms, problem, samp_opt);  % no. of start points given by samp_opt
        end
        samp_val_best=-samp_val_best;
    elseif strcmpi(samp_opt,'svm')
        % Quadratic SVM with class-balanced accuracy and value-based costs

        % Build training set
        X = [samp_1; samp_2];                        % rows = observations
        Y = [ones(size(samp_1,1),1);                 % +1 for samp_1
            -ones(size(samp_2,1),1)];               % -1 for samp_2

        % Class-balanced sample weights: each class gets total weight 0.5
        N1 = size(samp_1,1);
        N2 = size(samp_2,1);
        W  = [0.5/N1 * ones(N1,1);
            0.5/N2 * ones(N2,1)];

        % vals(i,j): value of assigning a point from sample i as j (i,j in {1,2})
        % Convert to misclassification costs (assuming vals(i,i) >= vals(i,≠i))
        c12 = vals(1,1) - vals(1,2);   % true class 1 → predicted 2
        c21 = vals(2,2) - vals(2,1);   % true class 2 → predicted 1

        % Cost matrix is ordered according to ClassNames = [-1 1]:
        %   row = true class, col = predicted class
        %   row/col 1: class -1 (sample 2)
        %   row/col 2: class +1 (sample 1)
        Cost = [0    c21;   % true -1, predicted +1
            c12  0  ];  % true +1, predicted -1

        % Train quadratic SVM
        svm_model = fitcsvm(X, Y, ...
            'KernelFunction',  'polynomial', ...
            'PolynomialOrder', 2, ...
            'Standardize',     true, ...
            'Weights',         W, ...
            'Cost',            Cost, ...
            'ClassNames',      [-1 1]);

        % Extract quad boundary coeffs

        % --- 1) Quadratic form in standardized+scaled space (u-space) ---
        alpha = svm_model.Alpha;                    % s x 1
        ysv   = svm_model.SupportVectorLabels;      % s x 1  (±1)
        A     = alpha .* ysv;               % A_i = alpha_i * y_i

        Z = svm_model.SupportVectors;               % standardized support vectors (z-space)
        s = svm_model.KernelParameters.Scale;       % kernel scale (scalar)

        U = Z ./ s;                         % u-space: u = z / s

        Q2_u = U' * (diag(A) * U);          % sum_i A_i u_i u_i^T
        q1_u = 2 * (U' * A);                % 2 * sum_i A_i u_i
        q0_u = sum(A) + svm_model.Bias;             % sum_i A_i + bias

        Q2_u = 0.5 * (Q2_u + Q2_u');        % enforce symmetry numerically

        % --- 2) Transform to raw x-space: u = D (x - mu) ---

        mu    = svm_model.Mu(:);                    % means used for standardization
        Sigma = svm_model.Sigma(:);                 % std devs used for standardization
        D     = diag(1 ./ (Sigma * s));     % u = D (x - mu)

        q2 = D' * Q2_u * D;

        samp_bd_svm.q2 = 0.5 * (q2 + q2');        % symmetrize again just in case
        samp_bd_svm.q1 = -2 * q2 * mu + D' * q1_u;
        samp_bd_svm.q0 = mu' * q2 * mu - q1_u' * D * mu + q0_u;

        samp_bd_current=[samp_bd_svm.q2(triu(true(size(samp_bd_svm.q2)))); samp_bd_svm.q1(:); samp_bd_svm.q0];


        % Compute class-balanced accuracy on training data
        Yhat     = predict(svm_model, X);
        hit_rate = mean(Yhat(Y ==  1) ==  1);   % TPR
        corr_rej = mean(Yhat(Y == -1) == -1);   % TNR
        samp_val_best = 0.5 * (hit_rate + corr_rej);

    elseif strcmpi(samp_opt,'smooth') % if smoothened accuracy, use fminunc
        if plots
            colors=colororder;
            hold on
            plot_sample(samp_1,.5,colors(1,:)); % HERE THE PRIORS ARE SET BY HAND TO 0.5, CHANGE THIS TO ACCEPT THE ACTUAL PRIOR
            plot_sample(samp_2,.5,colors(2,:));
            axis image; % axis([-10 10 -10 10])

            bd_handle=plot_boundary(norm_bd,dim,'plot_type','line','line_color',[0 1 0]); % placeholder for sample boundary
            plot_boundary(norm_bd,dim,'plot_type','line'); % plot the normal boundary
        end

        for k=[10.^linspace(-3,3,20) inf] % k is the sharpness of the accuracy function
            obj_fun=@(x) -samp_value_flat(x,samp_1,samp_2,'acc_sharpness',k,varargin{:}); % smooth objective function

            % optimize starting from current optimal samp_bd:
            samp_bd_opt_current=fminunc(obj_fun,samp_bd_current,optimset('Display','notify','TolX',0,'TolFun',1/(size(samp_1,1)+size(samp_2,1)),'MaxFunEvals',1e3*length(samp_bd_current)));

            % optimize starting from norm_bd:
            samp_bd_opt_norm=fminunc(obj_fun,norm_bd_flat,optimset('Display','notify','TolX',0,'TolFun',1/(size(samp_1,1)+size(samp_2,1)),'MaxFunEvals',1e3*length(norm_bd_flat)));

            % if this optimized boundary gives better exact classification,
            % check which one gives better performance:
            samp_val_opt_current=samp_value_flat(samp_bd_opt_current,samp_1,samp_2,varargin{:});
            samp_val_opt_norm=samp_value_flat(samp_bd_opt_norm,samp_1,samp_2,varargin{:});
            [samp_val_best,best_idx]=max([samp_val_current,samp_val_opt_current,samp_val_opt_norm]);

            if best_idx~=1 % if either of the two optimizations yielded a better result
                if best_idx==2 % if optimizing from current is the best
                    %             samp_val=samp_val_opt_current;
                    samp_bd_current=samp_bd_opt_current; % set it as the new optimal boundary
                elseif best_idx==3 % if optimizing from norm_bd is the best
                    samp_bd_current=samp_bd_opt_norm; % set it as the new optimal boundary
                end

                % make the boundary into a struct and plot
                q2=zeros(dim);
                q2(triu(true(dim)))=samp_bd_current(1:(dim^2+dim)/2);
                q2=q2+triu(q2,1)';
                samp_bd.q2=q2;
                samp_bd.q1=samp_bd_current(end-dim:end-1);
                samp_bd.q0=samp_bd_current(end);

                if plots
                    pause;
                    bd_handle.Function=quad2fun(samp_bd,1);
                    title(sprintf('sharpness = %g',k))
                end

                %                 fprintf('k=%.2f: %d %d %d %d \n',[k, samp_val_current,samp_val_opt_current,samp_val_opt_norm samp_val_best])
                samp_val_current=samp_val_best;
            end
        end
    end
    fprintf('Found boundary to optimize sample classification accuracy/value: \n %g (with normal boundary) \n %g (with sample boundary)\n',[samp_val_current, samp_val_best])
end