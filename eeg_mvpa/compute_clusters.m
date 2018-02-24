function [posSizes, negSizes, posLabels, negLabels] = compute_clusters(obsData,pVals,indiv_pval)

% step 1, initialize clusters
signCluster = obsData > 0;
signCluster(obsData < 0) = -1;

% step 2, apply threshold and determine positive and negative clusters
clusterMatrix = zeros(size(pVals));
clusterMatrix(pVals < indiv_pval) = 1;
posClusters = zeros(size(clusterMatrix));
negClusters = zeros(size(clusterMatrix));
posClusters(signCluster > 0) = clusterMatrix(signCluster > 0);
negClusters(signCluster < 0) = clusterMatrix(signCluster < 0);

% step 3, label clusters
posLabels = bwlabel(posClusters);
negLabels = bwlabel(negClusters);

% step 4 compute the sum of the perf stats in each of the clusters, separately
% for the positive and negative clusters
labels = 1:max(unique(posLabels));
for c = 1:numel(labels)
    posSizes(c) = sum(sum(obsData(posLabels==labels(c))));
end
labels = 1:max(unique(negLabels));
for c = 1:numel(labels)
    negSizes(c) = abs(sum(sum(obsData(negLabels==labels(c)))));
end
if ~exist('posSizes','var')
    posSizes = 0;
end
if ~exist('negSizes','var')
    negSizes = 0;
end