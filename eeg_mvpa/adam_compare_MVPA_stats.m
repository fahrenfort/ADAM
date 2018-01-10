function difstats = adam_compare_MVPA_stats(cfg,stats1,stats2,mask)
% ADAM_COMPARE_MVPA_STATS performs a statistical comparison between two stats structures that result
% from ADAM_COMPUTE_GROUP_MVPA or ADAM_COMPUTE_GROUP_ERP, given that the stats have the same
% dimensions (e.g. two time-time generalization matrices of equal size and same number of
% participants; i.e. within-subject repeated measures design).
%
% Use as:
%   adam_compare_MVPA_stats(cfg,stats1,stats2,mask);
%
% The function effectively does a condition subtraction (stats2 - stats1) and tests against zero
% using the same statistical procedure options as in adam_compute_group_MVPA. Each stats structure
% can be a 1xN structure array, where each corresponding element of the array is compared accross
% the two stat structures (see example below).
%
% Optionally, you can provide a mask: a binary matrix (for time-time or time-frequency) or vector
% (for ERP or MVPA with reduced_dim) to pre-select a 'region of interest' to constrain the
% comparison. 'mask' should be the fourth input argument. You can for example base the mask on a
% statistical outcome of the adam_compute_group_ functions or, extracted from the stats.pVals (see
% example below).
%
% The cfg (configuration) input structure can contain the following:
%
%       cfg.mpcompcor_method = 'uncorrected' (default); string specifying the method for multiple
%                              correction correction; other options are: 'cluster_based' for
%                              cluster-based permutation testing, 'fdr' for false-discovery rate,
%                              or 'none' if you don't wish to perform a statistical analysis.
%       cfg.indiv_pval       = .05 (default); integer; the statistical threshold for each individual
%                              time point; the fdr correction is applied on this threshold.
%       cfg.cluster_pval     = .05 (default); integer; if mpcompcor_method is set to
%                              'cluster_based', this is the statistical threshold for evaluating
%                              whether a cluster of significant contiguous time points (after the
%                              indiv_pval threshold) is larger than can be expected by chance; the
%                              cluster_pval should never be higher than the indiv_pval.
%       cfg.tail             = 'both' (default); string specifiying whether the statistical tests
%                              are done right- ('right') left- ('left'), or two-tailed ('both').
%                              Right-tailed tests for positive values, left-tailed tests for
%                              negative values, two-tailed tests test for both positive and negative
%                              values.
%
% The output diffstats structure will contain the following fields:
%
%       stats.ClassOverTime:        NxM matrix; group-average classification accuracy over N 
%                                   testing time points and M training time points; note that if
%                                   reduce_dims is specified, M will be 1, and ClassOverTime
%                                   will be squeezed to a Nx1 matrix of classification over time.
%                                   Here, ClassOverTime will be the difference in classification
%                                   accuracy between two conditions.
%       stats.StdError:             NxM matrix; standard-error across subjects over time-time
%       stats.pVals:                NxM matrix; p-values of each tested time-time point
%       stats.pStruct:              struct; cluster info, currently only appears if mpcompcor_method
%                                   was set to 'cluster_based'
%       stats.settings:             struct; the settings used during the level-1 (single subject)
%                                   results
%       stats.condname:             string; name of format 'condition1 - condition2'.
%       stats.cfg:                  struct; the cfg used to create these stats
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Example usage: 
%
% (1)
% mask = stats1.pVals'<1; 
%
% --> the pVals matrix consists of significant p-values within significant clusters, and 1's for all
%     non-signficant time(-time) points; so evaluating against <1 will give a binary mask with 1's
%     for significant, and 0's for non-significant pixels; this binary mask can be used for
%     comparing conditions:
%
% cfg = [];
% diffstats = adam_compare_MVPA_stats(cfg,stats1,stats2,mask);
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
% part of the ADAM toolbox, by J.J.Fahrenfort, VU, 2017/2018
% 
% See also ADAM_COMPUTE_GROUP_MVPA, ADAM_MVPA_FIRSTLEVEL, ADAM_PLOT_BDM_WEIGHTS

% input checking
if numel(stats1) ~= numel(stats2)
    error('The two stats input variables must contain the same number of elements.');
end
% get some defaults
if nargin<4
    mask = [];
end
reduce_dims = '';
tail = 'both';
cluster_pval = .05;
indiv_pval = .05;
mpcompcor_method = 'uncorrected';
% unpack original cfg
if isfield(stats1,'cfg')
    v2struct(stats1(1).cfg,{'fieldNames','reduce_dims','tail','cluster_pval','indiv_pval','mpcompcor_method','trainlim','testlim','reduce_dims'});
elseif isfiled(stats2,'cfg')
    v2struct(stats2(1).cfg,{'fieldNames','reduce_dims','tail','cluster_pval','indiv_pval','mpcompcor_method','trainlim','testlim','reduce_dims'});
end
v2struct(cfg);
if exist('one_two_tailed','var')
    error('The cfg.one_two_tailed field has been replaced by the cfg.tail field. Please replace cfg.one_two_tailed with cfg.tail using ''both'', ''left'' or ''right''. See help for further info.');
end

% pack cfg with defaults
nameOfStruct2Update = 'cfg';
cfg = v2struct(reduce_dims,tail,cluster_pval,indiv_pval,tail,mpcompcor_method,trainlim,testlim,reduce_dims,nameOfStruct2Update);

difstats = [];
for cStats = 1:numel(stats1)
    difstats = [difstats sub_compare_MVPA_stats(cfg,stats1(cStats),stats2(cStats),mask)];
end

function difstats = sub_compare_MVPA_stats(cfg,stats1,stats2,mask)
% unpack cfg
v2struct(cfg);
% compute some values
ClassTotal{1} = stats1.indivClassOverTime;
ClassTotal{2} = stats2.indivClassOverTime;
nSubj = size(ClassTotal{1},1);
%difstats.ClassOverTime = squeeze(mean(ClassTotal{1}-ClassTotal{2}))';
%difstats.StdError = squeeze(std(ClassTotal{1}-ClassTotal{2})/sqrt(nSubj))';
difstats.ClassOverTime = shiftdim(mean(ClassTotal{1}-ClassTotal{2}));
difstats.StdError = shiftdim(std(ClassTotal{1}-ClassTotal{2})/sqrt(nSubj));
difstats.condname = [stats1.condname ' - ' stats2.condname];
settings = stats1.settings; % assuming these are the same!
settings.measuremethod = 'accuracy difference';

% mask size
if isempty(mask)
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
difstats.pVals = ClassPvals;
difstats.pStruct = pStruct;
difstats.cfg = cfg;
difstats.settings = settings;