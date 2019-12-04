% A.Crimi 2019
% 
% Code related to Kara et al. 2020
%
function plot_heatmap(single_hm, feature_name)

% Plot heatmap 
heatmap = reshape(single_hm,24,16);
figure; imagesc(heatmap');
yticklabels={'B','D','F','H','J','L','N','P'};
set(gca,'ytick',[2     4     6     8    10    12    14    16],'yticklabel',yticklabels);
colorbar;
title(feature_name);
saveas( gcf, feature_name, 'png' );
