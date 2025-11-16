function [samp_val,samp_val_mat,samp_1_correct,samp_2_correct,samp_1_dv,samp_2_dv]=samp_value(samp_1,samp_2,dom,varargin)
% Expected value given two samples and a boundary. If outcome values
% are not additionally specified, this is the classification accuracy.
% Credits:
%   Abhranil Das <abhranil.das@utexas.edu>
%	Wilson S Geisler
%	Center for Perceptual Systems, University of Texas at Austin
% If you use this code, please cite:
%   A new method to compute classification error
%   https://jov.arvojournals.org/article.aspx?articleid=2750251

% parse inputs
parser=inputParser;
parser.KeepUnmatched=true;
addRequired(parser,'samp_1',@isnumeric);
addRequired(parser,'samp_2',@isnumeric);
addRequired(parser,'dom', @(x) isstruct(x) || isa(x,'function_handle'));
addParameter(parser,'dom_type','quad');
addParameter(parser,'d_scale_type','squeeze_dv', @(s) strcmpi(s,'squeeze_dv') || strcmpi(s,'squeeze_dist'));
addParameter(parser,'d_scale',1);
addParameter(parser,'acc_sharpness',inf,@isscalar); % sharpness of sigmoidal accuracy function. Inf means exact step function.
addParameter(parser,'vals',eye(2), @(x) isnumeric(x) && ismatrix(x));
addParameter(parser,'samp_balance',false,@islogical);

parse(parser,samp_1,samp_2,dom,varargin{:});
dom_type=parser.Results.dom_type;
acc_sharpness=parser.Results.acc_sharpness;
vals=parser.Results.vals;
d_scale_type=parser.Results.d_scale_type;
d_scale=parser.Results.d_scale;
samp_balance=parser.Results.samp_balance;

if strcmpi(dom_type,'ray_trace')
    [~,~,samp_1_correct]=dom(samp_1',[]);
    [~,~,samp_2_correct]=dom(samp_2',[]);
    samp_2_correct=~samp_2_correct;
else

    % compute sample decision variables
    if strcmpi(dom_type,'quad')
        q2=dom.q2;
        q1=dom.q1;
        q0=dom.q0;

        samp_1_dv=dot(samp_1,samp_1*q2',2) + samp_1*q1 + q0;
        samp_2_dv=dot(samp_2,samp_2*q2',2) + samp_2*q1 + q0;

    elseif strcmpi(dom_type,'fun')
        samp_1_cell=num2cell(samp_1,1);
        samp_1_dv=dom(samp_1_cell{:});
        samp_2_cell=num2cell(samp_2,1);
        samp_2_dv=dom(samp_2_cell{:});
    end

    % scale d' by squeezing decision variables together, so their
    % medians come to 0
    if strcmpi(d_scale_type,'squeeze_dv') && d_scale ~=1
        samp_1_dv=samp_1_dv-median(samp_1_dv)*(1-d_scale);
        samp_2_dv=samp_2_dv-median(samp_2_dv)*(1-d_scale);
    end

    % compute correct trials by binarizing dv
    if isinf(acc_sharpness)
        % step accuracy function:
        samp_1_correct=samp_1_dv > 0;
        samp_2_correct=samp_2_dv < 0;
    else
        % smoothened sigmoidal accuracy function:
        samp_1_correct=1./(1+exp(-acc_sharpness*samp_1_dv));
        samp_2_correct=1./(1+exp(acc_sharpness*samp_2_dv));
    end
end

samp_count_mat=[sum(samp_1_correct) sum(~samp_1_correct);
                sum(~samp_2_correct)  sum(samp_2_correct)];

samp_val_mat=samp_count_mat.*vals;

if ~samp_balance % if it is not required to be class-balanced
    samp_val=sum(samp_val_mat(:));
else
    % class-balanced expected value per sample point, i.e. average of the expected values per point in
    % each class
    samp_val=mean(sum(samp_val_mat,2)./sum(samp_count_mat,2));
end

end


