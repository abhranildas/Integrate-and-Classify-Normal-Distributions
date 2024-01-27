function samp_ov_err=histogram_overlap_err(samp_1,samp_2)

% compute the Bayes-optimal error using
% the overlap area of two histograms

dim=size(samp_1,2);
if dim==1
    [~,edges]=histcounts([samp_1;samp_2]);
    bin_width=diff(edges); bin_width=bin_width(1);

    h_a=histcounts(samp_1,edges,'normalization','pdf');
    h_b=histcounts(samp_2,edges,'normalization','pdf');

    samp_ov_err=sum(min(h_a,h_b))*bin_width/2;
elseif dim==2
    [~,x_edges,y_edges]=histcounts2([samp_1(:,1);samp_2(:,1)],[samp_1(:,2);samp_2(:,2)]);
    x_bin_width=diff(x_edges); x_bin_width=x_bin_width(1);
    y_bin_width=diff(y_edges); y_bin_width=y_bin_width(1);

    h_a=histcounts2(samp_1(:,1),samp_1(:,2),x_edges,y_edges,'normalization','pdf');
    h_b=histcounts2(samp_2(:,1),samp_2(:,2),x_edges,y_edges,'normalization','pdf');

    samp_ov_err=sum(min(h_a(:),h_b(:)))*x_bin_width*y_bin_width/2;
end

