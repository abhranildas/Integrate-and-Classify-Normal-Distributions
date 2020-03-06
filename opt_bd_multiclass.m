function [r,r_sign]=opt_bd_multiclass(n,normals,varargin)
% Return distances and signs (relative to mu_a) of the optimal quadratic
% boundary between normals a and b, in the unit direction n.
%
% How to use this command:
% See github readme at https://github.com/abhranildas/classify
%
% Credits:
%   Abhranil Das <abhranil.das@utexas.edu>
%	Wilson S Geisler
%	Center for Perceptual Systems, University of Texas at Austin
% If you use this code, please cite:
%   A new method to compute classification error
%   https://jov.arvojournals.org/article.aspx?articleid=2750251

% parse inputs
parser = inputParser;
addParameter(parser,'vals',[], @(x) isnumeric(x) && ismatrix(x)); % TODO default vals
parse(parser,varargin{:});

n_normals=length(normals);
n_dirs=size(n,2);
%vals=parser.Results.vals;

% WHAT TO RETURN AT ORIGIN

all_roots={};
all_signs={};
for i=2:n_normals
    % find boundary distance and sign
    [r,r_sign]=opt_bd(n,[normals(1).mu, normals(1).v],[normals(i).mu, normals(i).v],...
        'p_a',normals(1).p/(normals(1).p+normals(i).p));
    all_roots=vertcat(all_roots,r);
    all_signs=vertcat(all_signs,r_sign);
    %all_roots(i-1).r_sign=r_sign;
end

% merge roots from all the boundaries INEFFICIENT CODE
r={}; r_sign={};
for i=1:n_dirs
    r{i}=horzcat(all_roots{:,i});
    r_sign{i}=horzcat(all_signs{:,i});
end

% sort roots and signs
[r,idx]=cellfun(@sort,r,'un',0);
r_sign=cellfun(@(x,y) x(y),r_sign,idx,'un',0);

for i=1:n_dirs % INEFFICIENT CODE looping through the cell
    % make matrix of boundary signs at roots 
    signmat=-ones(n_normals-1,length(r{i}));
    for j=1:n_normals-1
        if ~isempty(all_roots{j,i}) % if this boundary exists in this direction
            for k=1:length(all_roots{j,i})
                % set root location to 0
                root=all_roots{j,i}(k);
                root_sign=all_signs{j,i}(k);
                root_idx=find(r{i}==root);
                signmat(j,root_idx)=0;
                
                % set signs of roots to left and right
                if sign(root*root_sign)>0
                    signmat(j,root_idx+1:end)=1;
                else
                    signmat(j,1:root_idx-1)=1;
                end
            end
        else % if this boundary doesn't exist in this direction
            signmat(j,:)=1;
        end
    end
    % reduce r and r_sign (to the ones corr. to the matrix columns >=0)
    r{i}=r{i}(all(signmat>=0));
    r_sign{i}=r_sign{i}(all(signmat>=0));
end