function stats = adam_compare_MVPA_stats(cfg,stats1,stats2,mask)
% function stats = compare_MVPA_stats(cfg,stats1,stats2,mask)
% compares and subtracts stats for plotting: stats1 - stats2
% cfg determines how statistics are computed
% cfg.tail = 'both' (can also be 'right' or 'left')
% cfg.indiv_pval = .05;
% cfg.cluster_pval = .05;
% cfg.mpcompcor_method = 'uncorrected' (default, can also be 'cluster_based', 'fdr' or 'none')

% get some defaults
reduce_dims = '';
tail = 'both';
cluster_pval = .05;
indiv_pval = .05;
mpcompcor_method = 'uncorrected';
% unpack original cfg
if isfield(stats1,'cfg')
    v2struct(stats1.cfg,{'fieldNames','reduce_dims','tail','cluster_pval','indiv_pval','mpcompcor_method','trainlim','testlim','reduce_dims'});
elseif isfiled(stats2,'cfg')
    v2struct(stats2.cfg,{'fieldNames','reduce_dims','tail','cluster_pval','indiv_pval','mpcompcor_method','trainlim','testlim','reduce_dims'});
end
v2struct(cfg);
if exist('one_two_tailed','var')
    error('The cfg.one_two_tailed field has been replaced by the cfg.tail field. Please replace cfg.one_two_tailed with cfg.tail using ''both'', ''left'' or ''right''. See help for further info.');
end

% pack cfg with defaults
nameOfStruct2Update = 'cfg';
cfg = v2struct(reduce_dims,tail,cluster_pval,indiv_pval,tail,mpcompcor_method,trainlim,testlim,reduce_dims,nameOfStruct2Update);

% compute some values
ClassTotal{1} = stats1.indivClassOverTime;
ClassTotal{2} = stats2.indivClassOverTime;
nSubj = size(ClassTotal{1},1);
stats.ClassOverTime = squeeze(mean(ClassTotal{1}-ClassTotal{2}))';
stats.StdError = squeeze(std(ClassTotal{1}-ClassTotal{2})/sqrt(nSubj))';
stats.condname = [stats1.condname ' - ' stats2.condname];
settings = stats1.settings; % assuming these are the same!
settings.measuremethod = 'accuracy difference';

% statistical testing
if nargin < 4
	mask = ones([size(ClassTotal{1},2) size(ClassTotal{1},3)]);
end
    
if strcmpi(mpcompcor_method,'fdr')
    % FDR CORRECTION
    [~,ClassPvals] = ttest(ClassTotal{1},ClassTotal{2},indiv_pval,'tail',tail);
    thresh = fdr(squeeze(ClassPvals),cluster_pval);
    ClassPvals(ClassPvals>thresh) = 1;
    pStruct = []; % still need to implement
elseif strcmpi(mpcompcor_method,'cluster_based')
    % CLUSTER BASED CORRECTION
    [ ClassPvals, pStruct ] = cluster_based_permutation(ClassTotal{1},ClassTotal{2},cfg,settings,mask);
elseif strcmpi(mpcompcor_method,'uncorrected')
    % NO MP CORRECTION
    [~,ClassPvals] = ttest(ClassTotal{1},ClassTotal{2},'tail',tail);
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
stats.cfg = cfg;
stats.settings = settings;