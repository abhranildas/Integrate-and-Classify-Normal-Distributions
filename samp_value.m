function [samp_val,samp_val_mat,samp_1_correct,samp_2_correct]=samp_value(samp_1,samp_2,dom,varargin)
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
parser=inputParser;
parser.KeepUnmatched=true;
addRequired(parser,'samp_1',@isnumeric);
addRequired(parser,'samp_2',@isnumeric);
addRequired(parser,'dom', @(x) isstruct(x) || isa(x,'function_handle'));
addParameter(parser,'dom_type','quad');
addParameter(parser,'vals',eye(2), @(x) isnumeric(x) && ismatrix(x));

parse(parser,samp_1,samp_2,dom,varargin{:});
dom_type=parser.Results.dom_type;
vals=parser.Results.vals;

if strcmpi(dom_type,'quad')
    q2=dom.q2;
    q1=dom.q1;
    q0=dom.q0;
    
    samp_1_correct=dot(samp_1,samp_1*q2',2) + samp_1*q1 + q0 > 0;
    samp_2_correct=dot(samp_2,samp_2*q2',2) + samp_2*q1 + q0 < 0;
    
elseif strcmpi(dom_type,'ray_scan')    
    [~,~,samp_1_correct]=dom(samp_1',[],[]);
    [~,~,samp_2_correct]=invert_dom(dom,samp_2');
    
elseif strcmpi(dom_type,'fun')
    samp_1_cell=num2cell(samp_1,1);
    samp_1_correct=dom(samp_1_cell{:})>0;
    samp_2_cell=num2cell(samp_2,1);
    samp_2_correct=dom(samp_2_cell{:})<0;
end

samp_val_mat=zeros(2);
samp_val_mat(1,1)=sum(samp_1_correct)*vals(1,1);
samp_val_mat(1,2)=sum(~samp_1_correct)*vals(1,2);
samp_val_mat(2,2)=sum(samp_2_correct)*vals(2,2);
samp_val_mat(2,1)=sum(~samp_2_correct)*vals(2,1);

samp_val=sum(samp_val_mat(:));

end


