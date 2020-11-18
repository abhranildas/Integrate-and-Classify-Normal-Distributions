function [merged_init_sign_ray,merged_x_ray]=combine_doms_ray(all_init_sign_ray,all_x_ray,op)
n_regs=length(all_init_sign_ray);

% merged initial sign
if strcmpi(op,'and')
merged_init_sign_ray=min(all_init_sign_ray); % -1 if any are -1, else 0 if any are 0, else 1 if all are 1.
elseif strcmpi(op,'or')
merged_init_sign_ray=max(all_init_sign_ray);
end

% merge and sort roots from all the boundaries
merged_x_ray=unique(horzcat(all_x_ray{:}));

% make matrix of boundary signs at roots
signmat=ones(n_regs,length(merged_x_ray));
for i=1:n_regs
    root_locs=ismember(merged_x_ray,all_x_ray{i});
    signmat(i,:)=(-1).^cumsum(root_locs);
    signmat(i,root_locs)=0; % set root locations to 0
end
signmat=signmat.*all_init_sign_ray;

% if a root is on multiple boundaries, pick it if all slope signs are the same
valid_multi_roots=true(1,length(merged_x_ray));
for j=find((sum(signmat==0)>1)&(~sum(signmat==-1))) % for every column j with no -1 and multiple 0s
    idx=find(~signmat(:,j));
    multi_root_signs=all_init_sign_ray(idx).*(-1).^sum(~signmat(idx,1:j),2);
    if ~(all(multi_root_signs==1) || all(multi_root_signs==-1)) % if signs are different
        valid_multi_roots(j)=0;
    end
end

if strcmpi(op,'and')
    % reduce x to the ones corr. to the sign matrix columns >=0, and valid multi roots
    merged_x_ray=merged_x_ray(all(signmat>=0) & valid_multi_roots);
elseif strcmpi(op,'or')
    % reduce x to the ones corr. to the sign matrix columns <=0, and valid multi roots
    merged_x_ray=merged_x_ray(all(signmat<=0) & valid_multi_roots);
end

if isequal(merged_x_ray,[])
    merged_x_ray=double.empty(1,0); % to allow matrix operations downstream
end