function difstats = adam_compare_MVPA_stats(cfg,stats1,stats2)
% ADAM_COMPARE_MVPA_STATS performs a statistical comparison between two stats structures that result
% from adam_compute_group_MVPA or adam_compute_group_ERP, given that the stats have the same
% dimensions (e.g. two time-time generalization matrices of equal size and same number of
% participants; i.e. within-subject repeated measures design).
%
% Use as:
%   adam_compare_MVPA_stats(cfg,stats1,stats2);
%
% The function effectively does a condition subtraction (stats2 - stats1) and tests against zero
% using the same statistical procedure options as in adam_compute_group_MVPA. Each stats structure
% can be a 1xN structure array, where each corresponding element of the array is compared accross
% the two stat structures (see example below).
%
% The cfg (configuration) input structure can contain the following:
%
%       cfg.mpcompcor_method = 'uncorrected' (default); string specifying the method for multiple
%                              correction correction; other options are: 'cluster_based' for
%                              cluster-based permutation testing, 'fdr' for false-discovery rate,
%                              or 'none' if you don't wish to perform a statistical analysis.
%       cfg.indiv_pval       = .05 (default); integer; the statistical threshold for each individual
%                              time point;
%       cfg.cluster_pval     = .05 (default); integer; if mpcompcor_method is set to
%                              'cluster_based', this is the statistical threshold for evaluating
%                              whether a cluster of contiguously significant time points (as
%                              determined by indiv_pval) is significant. If if mpcompcor_method is
%                              set to 'fdr', this is value that specifies the false discovery rate
%                              q (see help fdr_bh for details).
%       cfg.paired           = true; (default) for a paired t-test or false; for an unpaired
%                              t-test.
%       cfg.vartype          = 'equal' (default) whether an unpaired t-test assumes equal variances.
%       cfg.tail             = 'both' (default); string specifiying whether statistical tests are
%                              done right- ('right') or left-tailed ('left'), or two-tailed
%                              ('both'). Right-tailed tests for positive values, left-tailed tests
%                              for negative values, two-tailed tests test for both positive and
%                              negative values.
%       cfg.mask            =  Optionally, you can provide a mask: a binary matrix (for time-time
%                              or time-frequency) or vector (for ERP or MVPA with reduced_dim) to
%                              pre-select a 'region of interest' to constrain the comparison. You
%                              can for example base the mask on a statistical outcome of the
%                              adam_compute_group_ functions or, extracted from the stats.pVals
%                              (see example below).
%
% The output diffstats structure will contain the following fields:
%
%       stats.ClassOverTime:    NxM matrix; group-average classification accuracy over N
%                               testing time points and M training time points; note that if
%                               reduce_dims is specified, M will be 1, and ClassOverTime will be
%                               squeezed to a Nx1 matrix of classification over time. Here,
%                               ClassOverTime will be the difference in classification accuracy
%                               between two conditions.
%       stats.StdError:         NxM matrix; standard-error across subjects over time-time
%       stats.pVals:            NxM matrix; p-values of each tested time-time point 
%       stats.pStruct:          struct; cluster info, currently only appears if mpcompcor_method
%                               was set to 'cluster_based'
%       stats.settings:         struct; the settings used during the level-1 (single subject)
%                               results
%       stats.condname:         string; name of format 'condition1 - condition2'. 
%       stats.cfg:              struct; the cfg used to create these stats
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Example usage: 
%
% (1)
% cfg.mask = stats1.pVals<1; 
%
% --> the pVals matrix consists of significant p-values within significant clusters, and 1's for all
%     non-signficant time(-time) points; so evaluating against <1 will give a binary mask with 1's
%     for significant, and 0's for non-significant pixels; this binary mask can be used for
%     comparing conditions:
%
% cfg = [];
% diffstats = adam_compare_MVPA_stats(cfg,stats1,stats2);
%
% (2)
% cfg = [];
% diffstats = adam_compare_MVPA_stats(cfg,stats1,stats2);
%
% --> here, stats1 and stats2 could be EEG and MEG decoding, respectively, while each stats
%     structure has two class-comparisons (e.g. famous vs. non-famous faces and famous vs. scrambled
%     faces); diffstats will then also be a 1x2 structure-array reflecting the contrast EEG>MEG for:
%     famous vs non-famous, and famous vs scrambled.
%
% part of the ADAM toolbox, by J.J.Fahrenfort, VU/UvA, 2017/2018/2021
% 
% See also ADAM_COMPUTE_GROUP_MVPA, ADAM_MVPA_FIRSTLEVEL, ADAM_PLOT_MVPA, ADAM_PLOT_BDM_WEIGHTS, FDR_BH

% input checking
if numel(stats1) ~= numel(stats2)
    error('The two stats input variables must contain the same number of elements.');
end
% get some defaults
mask = [];
reduce_dims = '';
tail = 'both';
cluster_pval = .05;
indiv_pval = .05;
mpcompcor_method = 'uncorrected';
vartype = 'equal';
paired = true;
plot_dim = 'time_time';
trainlim = [];
testlim= [];
% unpack original cfg
if isfield(stats1,'cfg')
    v2struct(stats1(1).cfg);
elseif isfield(stats2,'cfg')
    v2struct(stats2(1).cfg);
end
v2struct(cfg);
if exist('one_two_tailed','var')
    error('The cfg.one_two_tailed field has been replaced by the cfg.tail field. Please replace cfg.one_two_tailed with cfg.tail using ''both'', ''left'' or ''right''. See help for further info.');
end

% pack cfg with defaults
nameOfStruct2Update = 'cfg';
cfg = v2struct(reduce_dims,tail,cluster_pval,indiv_pval,tail,mpcompcor_method,paired,vartype,trainlim,testlim,reduce_dims,plot_dim,nameOfStruct2Update);

difstats = [];
for cStats = 1:numel(stats1)
    difstats = [difstats sub_compare_MVPA_stats(cfg,stats1(cStats),stats2(cStats),mask)];
end

function difstats = sub_compare_MVPA_stats(cfg,stats1,stats2,mask)
% unpack cfg
vartype = 'equal';
v2struct(cfg);
% compute some values
ClassTotal{1} = stats1.indivClassOverTime;
ClassTotal{2} = stats2.indivClassOverTime;
% check whether to do paired or unpaired t-tests, in case someone did not specify
if paired && size(ClassTotal{2},1) ~= size(ClassTotal{1},1)
    paired = false;
    disp('not the same number of subjects, so doing an unpaired test');
end
if paired
    nSubj = size(ClassTotal{1},1);
    difstats.ClassOverTime = shiftdim(mean(ClassTotal{1}-ClassTotal{2}));
    difstats.indivClassOverTime = ClassTotal{1}-ClassTotal{2};
    difstats.StdError = shiftdim(std(ClassTotal{1}-ClassTotal{2})/sqrt(nSubj));
else
    nSubj1 = size(ClassTotal{1},1);
    nSubj2 = size(ClassTotal{2},1);
    difstats.ClassOverTime = shiftdim(mean(ClassTotal{1}) - mean(ClassTotal{2}));
    % equal variance
    difstats.StdError = shiftdim(sqrt((((nSubj1-1)*std(ClassTotal{1}).^2+(nSubj2-1)*std(ClassTotal{2}).^2)/(nSubj1+nSubj2-2)) * (1/nSubj1+1/nSubj2)));
    % unequal variance
    difstats.StdError = shiftdim(sqrt(std(ClassTotal{1}).^2/nSubj1+std(ClassTotal{2}).^2/nSubj2));
end
if ~isempty(strfind(stats1.condname,'subtract')) || ~isempty(strfind(stats1.condname,'average'))
    difstats.condname = [stats1.condname ' - ' stats2.condname];
else
    difstats.condname = ['subtract(' stats1.condname ',' stats2.condname ')'];
end
settings = stats1.settings; % assuming these are the same!
settings.measuremethod = [settings.measuremethod ' difference'];
settings.chance = 0;

% mask size
if isempty(mask)
	mask = ones([size(ClassTotal{1},2) size(ClassTotal{1},3)]);
end

% note that 'alpha' and 'tail' need to be specified explictly using name/value pairs, the ttest and
% ttest2 functions are no longer backwards compatible with Matlab 2012b
if any(strcmpi(mpcompcor_method,{'fdr','cluster_based','uncorrected'}))
    if paired
        [h,ClassPvals] = ttest(ClassTotal{1},ClassTotal{2},'alpha',indiv_pval,'tail',tail);
    else
        [h,ClassPvals] = ttest2(ClassTotal{1},ClassTotal{2},'alpha',indiv_pval,'tail',tail,'vartype',vartype);
    end
    ClassPvals = squeeze(ClassPvals);
else
    cfg.mpcompcor_method = 'none';
end
if strcmpi(mpcompcor_method,'fdr')
    % FDR CORRECTION
    pValsUncorrected = ClassPvals;
    h = fdr_bh(ClassPvals,cluster_pval,'dep');
    ClassPvals(~h) = 1;
    pStruct = compute_pstructs(h,ClassPvals,ClassTotal{1},ClassTotal{2},cfg,settings);
elseif strcmpi(mpcompcor_method,'cluster_based')
    % CLUSTER BASED CORRECTION
    pValsUncorrected = ClassPvals;
    [ ClassPvals, pStruct ] = cluster_based_permutation(ClassTotal{1},ClassTotal{2},cfg,settings,mask);
elseif strcmpi(mpcompcor_method,'uncorrected')
    % NO MCP CORRECTION
    pStruct = compute_pstructs(squeeze(h),ClassPvals,ClassTotal{1},ClassTotal{2},cfg,settings);
else
    % NO TESTING, PLOT ALL
    ClassPvals = zeros([size(ClassTotal{1},2) size(ClassTotal{1},3)]);
    pStruct = [];
end
cfg.chance = 0;
% output stats
difstats.pVals = ClassPvals;
if exist('pValsUncorrected','var')
    difstats.pValsUncorrected   = pValsUncorrected;
end
difstats.pStruct = pStruct;
difstats.cfg = cfg;
difstats.settings = settings;
% compute latency
try
    difstats.latencies = extract_latency(cfg,difstats);
catch ME
    disp('Cannot extract latencies.');
    disp(ME.message);
    difstats.latencies = [];
end