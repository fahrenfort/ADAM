function pStruct = compute_pstructs(H,pVals,data1,data2_or_chance_level,cfg,settings,mask)
% function computes pStruct start and stop time and/or frequency for both positive and negative
% clusters, and returns the result.

% compute mask
if nargin < 7 || isempty(mask)
    mask = ones([size(data1,2) size(data1,3)]);
elseif sum([size(data1,2) size(data1,3)] == [size(mask,1) size(mask,2)]) ~= ndims(mask) % ndims(mask) ~= ndims(data1)-1 || 
    error('mask does not have the same size as data1');
end
mask = logical(mask);

% repeat data if data2 is a chance variable
if numel(data2_or_chance_level) == 1
    data2 = repmat(data2_or_chance_level,size(data1));
else
    data2 = data2_or_chance_level;
end

% compute data sign
signCluster  = squeeze(mean(data1,1) - mean(data2,1));
  
% use mask to restrict relevant info
pVals(~mask) = 1;
signCluster(~mask) = 0;

% step 2, apply threshold and determine positive and negative clusters
posClusters = H;
negClusters = H;
posClusters(signCluster < 0) = 0;
negClusters(signCluster > 0) = 0;

% step 3, label clusters
posLabels = bwlabel(posClusters);
negLabels = bwlabel(negClusters);

% step 4, compute average p-values in pos and neg clusters
clusterPvals = pVals;
posVals = setdiff(unique(posLabels),0);
for cClust = 1:numel(posVals)
    clusterPvals(posLabels==posVals(cClust)) = mean(pVals(posLabels==posVals(cClust)));
end
negVals = setdiff(unique(negLabels),0);
for cClust = 1:numel(negVals)
    clusterPvals(negLabels==negVals(cClust)) = mean(pVals(negLabels==negVals(cClust)));
end

% compute pstructs of clusters
pStruct.posclusters = compute_pstruct(posLabels,clusterPvals,signCluster,cfg,settings);
pStruct.negclusters = compute_pstruct(negLabels,clusterPvals,-signCluster,cfg,settings);
