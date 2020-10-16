function plot_sample(dist,prior,color)
if ~exist('color','var')
    color='blue';
end
dim=size(dist,2);
if dim==1
    % plot sample histogram
    [heights,edges]=histcounts(dist,'BinMethod','scott','normalization','pdf');
    histogram('BinCounts',heights*prior,'BinEdges',edges,'facecolor',color,'facealpha',0.3,'edgecolor',color,'edgealpha',0.3);
elseif dim==2
    % plot sample points
    h=plot(dist(:,1),dist(:,2),'.','markersize',4,'color',hsv2rgb(rgb2hsv(color)-[0 .5 0]));
    uistack(h,'bottom')
    uistack(h,'up',3)
elseif dim==3
    % plot sample points
    plot3(dist(:,1),dist(:,2),dist(:,3),'.','color',color,'markersize',4);
end