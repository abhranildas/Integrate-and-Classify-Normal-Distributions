function [merged_init_sign_ray,merged_x_ray]=opt_reg_multi_ray(all_init_sign_ray,all_x_ray,n_normals)
% merged initial sign
merged_init_sign_ray=min(all_init_sign_ray); % -1 if any are -1, else 0 if any are 0, else 1 if all are 1.

% merge and sort roots from all the boundaries
merged_x_ray=unique(horzcat(all_x_ray{:}));

% make matrix of boundary signs at roots
signmat=ones(n_normals-1,length(merged_x_ray));
for i=1:n_normals-1
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

% reduce z to the ones corr. to the sign matrix columns >=0, and valid multi roots
merged_x_ray=merged_x_ray(all(signmat>=0) & valid_multi_roots);

if isequal(merged_x_ray,[])
    merged_x_ray=double.empty(1,0); % to allow matrix operations downstream
end