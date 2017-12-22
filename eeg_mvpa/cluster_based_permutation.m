function [ clusterPvals, pStruct ] = cluster_based_permutation(data1,data2_or_chance_level,cfg,settings,mask,connectivity)
% implements Maris, E., & Oostenveld, R. (2007). Nonparametric statistical
% testing of EEG- and MEG-data. Journal of Neuroscience Methods, 164(1),
% 177?190. http://doi.org/10.1016/J.Jneumeth.2007.03.024
%
% Note potential issues when using t-testing using accuracy measures to
% obtain population level inferences, e.g. read:
% Allefeld C, G?rgen K, Haynes J-D (2016) Valid population inference for
% information-based imaging: From the second-level t-test to prevalence
% inference. Neuroimage 141:378?392.
%
% data1 is a subject x dim1 x dim2 array, dim1 and dim2 can be any type of
% dimension (time, frequency, electrode etc) containing some dependent
% measure (e.g. classification accuracy, power values, anything).
% data2_or_chance either contains a data matrix of the same size as data1,
% or contains a single value corresponding to the value against which data1
% should be tested.
% cfg.indiv_pval and cfg.cluster_pval respectively determine
% the p-value used for inclusion into the cluser, and the p-value used for
% evaluation overall cluster significance
% cfg.paired (true or false) and cfg.tail ('one' or 'two)
% determine whether to apply paired t-testing and whether to apply one- or
% two-tailed testing
% cfg.iterations is the number of iterations in the permutation 
% mask is the same size as data1 without the subjects (so dim1 x dim2), and
% allows you to restrict the cluster based test to a specific region by
% specifying 0 for the regions that should not be included and 1 for the
% regions that should
% connectivity is a matrix that outlines which dim1 points are connected to
% which dim1 points (so size dim1 x dim1). Useful when doing a cluster
% based permutation test across electrodes.
%
% outputs clusterPvals with p-values < pval(2) for the significant
% clusters and 1's for all other points and outputs pStruct that shows
% where the clusters start and stop (either in the time domain and/or in
% the frequency domain)
%
% By J.J.Fahrenfort, VU, 2015, 2016

if nargin < 6
    connectivity = [];
end
if nargin < 5 || isempty(mask)
    mask = ones([size(data1,2) size(data1,3)]);
    disp(['using all ' num2str(numel(find(mask))) ' datapoints in cluster based permutation']);
elseif sum([size(data1,2) size(data1,3)] == [size(mask,1) size(mask,2)]) ~= ndims(mask) % ndims(mask) ~= ndims(data1)-1 || 
    error('mask does not have the same size as data1');
else
    disp(['there are ' num2str(numel(find(mask))) ' out of ' num2str(numel(mask)) ' data points in your mask during cluster based permutation']);
end
mask = logical(mask);

% some default settings (these are defined by cfg)
indiv_pval = .05; % default = .05
cluster_pval = .05; % default = .05
iterations = 1000;
paired = true;
tail = 'two';
v2struct(cfg);
tail = one_two_tailed; % jvd: added this line, otherwise it always did the default of two-tailed, also if one-tailed specified in cfg (22-12-17)
pval = [indiv_pval, cluster_pval];

% repeat data if data2 is a chance variable
if numel(data2_or_chance_level) == 1
    data2 = repmat(data2_or_chance_level,size(data1));
else
    data2 = data2_or_chance_level;
end

% check whether to do paired or unpaired t-tests
if size(data2,1) ~= size(data1,1)
    paired = false;
    disp('not the same number of subjects, so doing an unpaired test');
end

% step 1 to 5 compute observed cluster statistics
[actPosSizes, actNegSizes, posLabels, negLabels, sigClusters] = compute_cluster_sizes(data1,data2,pval(1),tail,mask,connectivity,paired);
clusterPvals = ones(size(sigClusters));

% step 6 iterate to determine how often permuted clusters exceed the observed cluster threshold
cPosClust = zeros(1,max(unique(posLabels)));
cNegClust = zeros(1,max(unique(negLabels)));

for cIt = 1:iterations
    if ~mod(cIt,round(iterations/10))
        fprintf(1,'iteration %d of %d\n', cIt, iterations);
    end
    
    % make random partitions
    if paired % for paired, keep observations paired under permutation
        randlabels = rand(1,size(data1,1))<.5;
        randdata1(randlabels,:,:) = data1(randlabels,:,:);
        randdata1(~randlabels,:,:) = data2(~randlabels,:,:);
        randdata2(randlabels,:,:) = data2(randlabels,:,:);
        randdata2(~randlabels,:,:) = data1(~randlabels,:,:);
    else % for unpaired, fully randomize observations under permutation
        alldata = cat(1,data1,data2);
        alldata = alldata(randperm(size(alldata,1)),:,:);
        randdata1 = alldata(1:size(data1,1),:,:);
        randdata2 = alldata(size(data1,1)+1:end,:,:);
    end

    % repeat step 1 to 5 recompute cluster sizes under random permutation
    [randPosSizes, randNegSizes] = compute_cluster_sizes(randdata1,randdata2,pval(1),tail,mask,connectivity,paired);
    maxRandSize = max([randPosSizes randNegSizes]);

    % count cluster p-values
    cPosClust = cPosClust + (maxRandSize > actPosSizes);
    cNegClust = cNegClust + (maxRandSize > actNegSizes);
end

% compute pvalues for clusters
pPos = cPosClust / iterations;
pNeg = cNegClust / iterations;

% cut out clusters that do not cross the threshold
if strcmp(tail,'two')
    for c=1:max(unique(posLabels))
        if pPos(c) < pval(2)/2
            clusterPvals(posLabels == c) = pPos(c);
        else
            posLabels(posLabels == c) = 0;
        end
    end
    for c=1:max(unique(negLabels))
        if pNeg(c) < pval(2)/2
            clusterPvals(negLabels == c) = pNeg(c);
        else
            negLabels(negLabels == c) = 0;
        end
    end
else
    for c=1:max(unique(posLabels))
        if pPos(c) < pval(2)
            clusterPvals(posLabels == c) = pPos(c);
        else
            posLabels(posLabels == c) = 0;
        end
    end
end
pStruct.posclusters = compute_pstruct(posLabels,clusterPvals,squeeze(mean(data1)-mean(data2)),cfg,settings,mask,connectivity);
pStruct.negclusters = compute_pstruct(negLabels,clusterPvals,squeeze(mean(data2)-mean(data1)),cfg,settings,mask,connectivity);

function [posSizes, negSizes, posLabels, negLabels, pVals] = compute_cluster_sizes(data1,data2,indiv_pval,one_two_tailed,mask,connectivity,paired)
% step 1, determine 'actual' p values
if strcmpi(one_two_tailed,'two')
    tail = 'both';
else
    tail = 'right';
end
% restrict data
for c = 1:size(data1)
    maskdata1(c,:,:) = data1(c,mask);
end
for c = 1:size(data2)
    maskdata2(c,:,:) = data2(c,mask);
end

pVals = ones(size(mask));
tVals = zeros(size(mask));
if paired
    [~,pVals(mask),~,stats] = ttest(maskdata1,maskdata2,indiv_pval,tail);
else
    [~,pVals(mask),~,stats] = ttest2(maskdata1,maskdata2,indiv_pval,tail);    
end
tVals(mask) = squeeze(stats.tstat);

% initialize clusters
signCluster  = squeeze(mean(data1,1) - mean(data2,1));

% use mask to restrict relevant info
pVals(~mask) = 1;
signCluster(~mask) = 0;

% step 2, apply threshold and determine positive and negative clusters
clusterMatrix = squeeze(pVals < indiv_pval);
posClusters = zeros(size(clusterMatrix));
negClusters = zeros(size(clusterMatrix));
posClusters(signCluster > 0) = clusterMatrix(signCluster > 0);
negClusters(signCluster < 0) = clusterMatrix(signCluster < 0);

% step 3, label clusters
if isempty(connectivity)
    posLabels = bwlabel(posClusters);
    negLabels = bwlabel(negClusters);
else
    % slightly more complex to find clusters in topomap data
    elecs2do = posClusters;
    posLabels = zeros(size(clusterMatrix)); % 
    cClust = 0;
    for c=1:numel(elecs2do)
        if elecs2do(c) == 1 % only look for a new cluster if this electrode has not been looked at yet
            [clustlabels, elecs2do] = find_elec_clusters(elecs2do,posClusters,c,connectivity);
            if sum(clustlabels) > 1
                cClust = cClust + 1;
                posLabels(clustlabels) = cClust;
            end
        end
    end
    elecs2do = negClusters;
    negLabels = zeros(size(clusterMatrix)); % 
    cClust = 0;
    for c=1:numel(elecs2do)
        if elecs2do(c) == 1 % only look for a new cluster if this electrode has not been looked at yet
            [clustlabels, elecs2do] = find_elec_clusters(elecs2do,negClusters,c,connectivity);
            if sum(clustlabels) > 1
                cClust = cClust + 1;
                negLabels(clustlabels) = cClust;
            end
        end
    end
end

% step 4 compute the sum of the t-stats in each of the clusters, separately
% for the positive and negative clusters
labels = 1:max(unique(posLabels));
for c = 1:numel(labels)
    posSizes(c) = sum(sum(tVals(posLabels==labels(c))));
end
labels = 1:max(unique(negLabels));
for c = 1:numel(labels)
    negSizes(c) = abs(sum(sum(tVals(negLabels==labels(c)))));
end
if ~exist('posSizes','var')
    posSizes = 0;
end
if ~exist('negSizes','var')
    negSizes = 0;
end

function [clustlabels, elecs2do] = find_elec_clusters(elecs2do,sigelectrodes,cElec,connectivity)
% this function calls itself recursively to find all connected electrodes
% elecs2do: electrodes that have not yet been inspected
% sigelectrodes: electrodes that were initially signficant
% cElec: current seed electrode
% connectivity: of electrodes on cap, computed by get_connected_electrodes.m
connectedlabels = connectivity(cElec,:)'; % find electrodes connected to this electrode
clustlabels = connectedlabels & sigelectrodes & elecs2do; % find all significant electrodes connected to this electrode that have not yet been inspected
labelindex = find(clustlabels); % list their indices
for c=1:numel(labelindex) % run through these electrodes
    elecs2do(labelindex(c)) = 0; % look only once, set to zero after inspection
    % calls itself recusively to grow the cluster:
    [clustlabels2, elecs2do] = find_elec_clusters(elecs2do,sigelectrodes,labelindex(c),connectivity);
    clustlabels = clustlabels | clustlabels2; % increment
end
