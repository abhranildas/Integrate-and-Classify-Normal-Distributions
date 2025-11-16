function samp_bd_current=optimize_samp_bd(samp_1,samp_2,norm_bd,varargin)
    parser=inputParser;
    parser.KeepUnmatched=true;
    addRequired(parser,'samp_1',@isnumeric);
    addRequired(parser,'samp_2',@isnumeric);
    addRequired(parser,'norm_bd',@isstruct);
    addParameter(parser,'samp_opt',true, @(x) strcmpi(x,'step') || strcmpi(x,'smooth') || islogical(x));
    addParameter(parser,'samp_opt_plot','step', @(x) strcmpi(x,'step') || strcmpi(x,'smooth') || x==false);
    addParameter(parser,'plots',false, @islogical);

    parse(parser,samp_1,samp_2,norm_bd,varargin{:});
    samp_opt=parser.Results.samp_opt;
    plots=parser.Results.plots;

    dim=size(samp_1,2);

    % start from the normal boundary
    norm_bd_flat=[norm_bd.q2(triu(true(size(norm_bd.q2)))); norm_bd.q1(:); norm_bd.q0];
    samp_bd_current=norm_bd_flat;
    samp_val_current=samp_value_flat(samp_bd_current,samp_1,samp_2,varargin{:});

    if samp_opt % if exact accuracy, use fminsearch
        obj_fun=@(x) -samp_value_flat(x,samp_1,samp_2,varargin{:}); % step objective function
        [samp_bd_current,samp_val_best]=fminsearch(obj_fun,samp_bd_current,optimset('Display','notify','TolX',0,'TolFun',1/(size(samp_1,1)+size(samp_2,1))));
        fprintf('Sample classification accuracy/value optimized from %g to %g \n',[samp_val_current, -samp_val_best])

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