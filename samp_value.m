function [samp_ex_val,samp_val_mat,samp_1_correct,samp_2_correct]=samp_value(samp_1,samp_2,reg,varargin)
% Expected value given two samples and a boundary. If outcome values
% are not additionally specified, this is the classification accuracy.
% Credits:
%   Abhranil Das <abhranil.das@utexas.edu>
%	R Calen Walshe
%	Wilson S Geisler
%	Center for Perceptual Systems, University of Texas at Austin
% If you use this code, please cite:
%   A new method to compute classification error
%   https://jov.arvojournals.org/article.aspx?articleid=2750251

% parse inputs
parser = inputParser;
addRequired(parser,'samp_1',@isnumeric);
addRequired(parser,'samp_2',@isnumeric);
addRequired(parser,'reg', @(x) isstruct(x) || isa(x,'function_handle'));
addParameter(parser,'reg_type','quad');
addParameter(parser,'vals',eye(2), @(x) isnumeric(x) && ismatrix(x));

parse(parser,samp_1,samp_2,reg,varargin{:});
reg_type=parser.Results.reg_type;
vals=parser.Results.vals;

if strcmpi(reg_type,'quad')
    q2=reg.q2;
    q1=reg.q1;
    q0=reg.q0;
    
    samp_1_correct=dot(samp_1,samp_1*q2',2) + samp_1*q1 + q0 > 0;
    samp_2_correct=dot(samp_2,samp_2*q2',2) + samp_2*q1 + q0 < 0;
    
elseif strcmpi(reg_type,'ray_scan')    
    [~,~,samp_1_correct]=reg(samp_1',[],[]);
    %reg_inv=@(n) invert_reg(reg,n);
    [~,~,samp_2_correct]=invert_reg(reg,samp_2'); %reg_inv(samp_2');
    
elseif strcmpi(reg_type,'fun')
    samp_1_cell=num2cell(samp_1,1);
    samp_1_correct=reg(samp_1_cell{:})>0;
    samp_2_cell=num2cell(samp_2,1);
    samp_2_correct=reg(samp_2_cell{:})<0;
end

samp_val_mat=zeros(2);
samp_val_mat(1,1)=sum(samp_1_correct)*vals(1,1);
samp_val_mat(1,2)=sum(~samp_1_correct)*vals(1,2);
samp_val_mat(2,2)=sum(samp_2_correct)*vals(2,2);
samp_val_mat(2,1)=sum(~samp_2_correct)*vals(2,1);

samp_ex_val=sum(samp_val_mat(:))/(size(samp_1,1)+size(samp_2,1));

end


