function [ blurdum, map ] = showstatsTFR(dum, prob, ZLIM)
% showstatsTFR Take a stat map from ft and blur out nonsig parts of the TFR
map = [];
% right color map
cmap  = brewermap([],'RdBu');
colormap(cmap(end:-1:1,:));
% load('colormap170613.mat');
nsalpha=0.2;
pval = squeeze(prob);
thrcluster = zeros(size(pval));
thrcluster = thrcluster+nsalpha; % set alpha level for non-significant stuff
thrcluster(pval<0.05) = 1; % set alpha level for significant stuff
tmpcdat = (dum + -ZLIM(1)) * (size(colormap,1) / (-ZLIM(1) + ZLIM(2))); % create indices within clim range for colormap
% ind->rgb->hsv
rgbcdat = ind2rgb(uint8(floor(tmpcdat)), colormap);
hsvcdat = rgb2hsv(rgbcdat);
% change saturation values
hsvcdat(:,:,2) = hsvcdat(:,:,2) .* thrcluster;
% hsv->rgb ->  plot
blurdum = hsv2rgb(hsvcdat);
%[blurdum, map] = rgb2ind(blurdum,colormap);