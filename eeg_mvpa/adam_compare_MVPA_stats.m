function stats = adam_compare_MVPA_stats(stats1,stats2,cfg,mask)
% function stats = compare_MVPA_stats(stats1,stats2,cfg,mask)
% compares and subtracts stats for plotting: stats1 - stats2
% cfg determines how statistics are computed
% cfg.one_two_tailed = 'two' (can also be 'one')
% cfg.indiv_pval = .05;
% cfg.cluster_pval = .05;
% cfg.mpcompcor_method = 'uncorrected' (default, can also be 'cluster_based', 'fdr' or 'none')

one_two_tailed = 'two';
cluster_pval = .05;
indiv_pval = .05;
mpcompcor_method = 'uncorrected';
v2struct(cfg);
pval(1) = indiv_pval;
pval(2) = cluster_pval;

% compute some values
ClassTotal{1} = stats1.indivClassOverTime;
ClassTotal{2} = stats2.indivClassOverTime;
nSubj = size(ClassTotal{1},1);
stats.ClassOverTime = squeeze(mean(ClassTotal{1}-ClassTotal{2}))';
stats.StdError = squeeze(std(ClassTotal{1}-ClassTotal{2})/sqrt(nSubj))';
stats.condname = [stats1.condname '-' stats2.condname];
settings = stats1.settings; % assuming these are the same!

% statistical testing
if nargin < 4
	mask = ones([size(ClassTotal{1},2) size(ClassTotal{1},3)]);
end
    
if strcmp(mpcompcor_method,'fdr')
    % FDR CORRECTION
    if strcmp(one_two_tailed,'two')
        [~,ClassPvals] = ttest(ClassTotal{1},ClassTotal{2},'tail','both');
    else
        [~,ClassPvals] = ttest(ClassTotal{1},ClassTotal{2},'tail','right');
    end
    ClassPvals = squeeze(ClassPvals);
    h = fdr_bh(ClassPvals,pval(2),'pdep');
    ClassPvals(~h) = 1;
    pStruct = []; % still need to implement
elseif strcmp(mpcompcor_method,'cluster_based')
    % CLUSTER BASED CORRECTION
    [ ClassPvals, pStruct ] = cluster_based_permutation(ClassTotal{1},ClassTotal{2},cfg,settings,mask);
elseif strcmp(mpcompcor_method,'uncorrected')
    % NO MP CORRECTION
    if strcmp(one_two_tailed,'two')
        [~,ClassPvals] = ttest(ClassTotal{1},ClassTotal{2},'tail','both');
    else
        [~,ClassPvals] = ttest(ClassTotal{1},ClassTotal{2},'tail','right');
    end
    ClassPvals = squeeze(ClassPvals);
    ClassPvals(~mask) = 1;
    pStruct = []; % still need to implement
else
    % NO TESTING, PLOT ALL
    ClassPvals = zeros([size(ClassTotal{1},2) size(ClassTotal{1},3)]);
    pStruct = [];
end

% output stats
stats.pVals = ClassPvals;
stats.pStruct = pStruct;
stats.settings = settings;