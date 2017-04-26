function [weightsubtraction, pStruct] = compare_weights(subjweights1,subjweights2,labels1,labels2,gsettings)
% function weightsubtraction = compare_weights(weights,subjweights,gsettings)
% THIS FUNCTION NEEDS SOME WORK, NEEDS TO BE MADE MORE GENERIC AND IN LINE
% WITH OTHER INPUT-OUTPUT FUNCTIONS
% compares and subtracts weights for plotting: subjweights1 - subjweights2
% subjweights are obtained from the output plot_FEM_weights. They are
% subsequently extracted from indivWeights in that output, and have the
% dimensions (subjects x electrodes) chanlocs contains the labels
% associated with subjweights1 and subjweights2.
% gsettings.paired = true or false determines whether to do a paired or
% unpaired t-test.
% output weights can be plotted using plot_MVPA_weights
%
% JJF,VU, 2016

if nargin < 2
    help compute_weight_stats;
    error('needs two arguments');
end
% set some defaults
one_two_tailed = 'two';
paired = true;
indiv_pval = .05;
cluster_pval = .05;
iterations = 1000;
mpcompcor_method = 'uncorrected';
v2struct(gsettings);
pval(1) = indiv_pval;
pval(2) = cluster_pval;

if size(subjweights1,1) ~= size(subjweights1,2)
    paired = false;
    disp('unequal number of subjects, doing unpaired comparison');
end

% equalize dimensions and put electrodes in the same order
if numel(labels1) > numel(labels2)
    [~, chanindex] = intersect(labels1,labels2,'stable');
    labels = labels1(chanindex);
    subjweights1 = subjweights1(:,chanindex);
    disp('removing electrodes from subjweights1 to maintain the same size as subjweights2');
else
    [~, chanindex] = intersect(labels2,labels1,'stable');
    labels = labels2(chanindex);
    subjweights2 = subjweights2(:,chanindex);
    disp('removing electrodes from subjweights2 to maintain the same size as subjweights1');
end

% do some statistics
if strcmpi(mpcompcor_method, 'cluster_based')
    connectivity = get_connected_electrodes(labels);
    settings.channels = labels';
    [clusterPvals, pStruct] = cluster_based_permutation(subjweights1,subjweights2,gsettings,settings,[],connectivity);
elseif strcmpi(mpcompcor_method, 'uncorrected')
    if strcmpi(one_two_tailed,'two')
        tail = 'both';
    else
        tail = 'right';
    end
    if paired
        [~,clusterPvals] = ttest(subjweights1,subjweights2,'tail',tail);
    else
        [~,clusterPvals] = ttest2(subjweights1,subjweights2,'tail',tail);
    end
elseif strcmpi(mpcompcor_method, 'fdr')
    error('fdr not yet implemented yet');
    clusterPvals = [];
else
    clusterPvals = [];
end

plot_data = mean(subjweights1,1) - mean(subjweights2,1);
elecs = find(clusterPvals<cluster_pval);

[ pStruct.spearman.corr, pStruct.spearman.pval ] = corr(mean(subjweights1,1)',mean(subjweights2,1)','type','Spearman','tail','both');
[ pStruct.pearson.corr, pStruct.pearson.pval ] = corr(mean(subjweights1,1)',mean(subjweights2,1)','type','Pearson','tail','both');

chanlocdata = readlocs('plotting_1005.sfp','importmode','native');
[~, chanindex] = intersect({chanlocdata(:).labels},labels,'stable');
chanlocs = chanlocdata(chanindex);

% plotting weights
fh = figure;
set(fh, 'Position', get(0,'Screensize'));
set(gcf,'color','w');
title('condition 1 - condition 2','FontSize', 32);
subplot(2,4,1);
if strcmp(imgtype,'vec')
    topoplot_jjf(plot_data,chanlocs','maplimits',weightlim,'style','blank','electrodes','on','plotrad',.6,'nosedir','-Y','emarker2',{elecs,'o','k',10,1}); %
elseif strcmp(imgtype,'png')
    topoplot_jjf(plot_data,chanlocs','maplimits',weightlim,'style','map','electrodes','off','plotrad',.6,'hcolor','none','shading','interp','nosedir','-Y'); %
else
    topoplot_jjf(plot_data,chanlocs','maplimits',weightlim,'style','map','electrodes','ptslabels','plotrad',.6,'nosedir','-Y','emarker2',{elecs,'o','k',10,1}); 
    cbar('vert');
end

% return some stuff
weightsubtraction.labels = labels;
weightsubtraction.clusterPvals = clusterPvals;
weightsubtraction.avWeights = plot_data;
