function samp_ov_err=histogram_overlap_err(samp_1,samp_2)

    % compute the Bayes-optimal error using
    % the overlap area of two histograms

    [~,edges]=histcounts([samp_1;samp_2]);
    bin_width=diff(edges); bin_width=bin_width(1);

    h_a=histcounts(samp_1,edges,'normalization','pdf');
    h_b=histcounts(samp_2,edges,'normalization','pdf');
    samp_ov_err=sum(min(h_a,h_b))*bin_width/2;