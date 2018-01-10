function stats = adam_compare_MVPA_stats(cfg,stats1,stats2,mask)
% ADAM_COMPARE_MVPA_STATS performs a statistical comparison between two stats structures that result
% from ADAM_COMPUTE_GROUP_MVPA, given that the stats have the same dimensions (e.g. two time-time
% generalization matrices of equal size and same number of participants; i.e. within-subject
% repeated measures design).
%
% Use as:
%   adam_compare_MVPA_stats(stats1,stats2,cfg,varargin);
%
% The function effectively does a condition subtraction (stats2 - stats1) and tests against zero
% using the same statistical procedure options as in ADAM_COMPUTE_GROUP_MVPA. 
%
% Optionally, you can provide a mask: a binary matrix (for time-time or time-frequency) or vector
% (for ERP or MVPA with reduced_dim) to pre-select a 'region of interest' to constrain the
% comparison. 'mask' should be the fourth input argument. You can for example base the mask on a
% statistical outcome of ADAM_COMPUTE_GROUP_MVPA, extracted from the stats.pVals (see example
% below).
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
%       cfg.tail             = 'both' (default); string specifiying whether the t-tests are done
%                              right- ('right') or left-tailed ('left'), or two-tailed ('both').
%
% The output diffstats structure will contain the following fields:
%
%       stats.ClassOverTime:        NxM matrix; group-average classification accuracy over N 
%                                   testing time points and M training time points; note that if
%                                   reduce_dims is specified, M will be 1, and ClassOverTime
%                                   will be squeezed to a Nx1 matrix of classification over time.
%                                   Here, ClassOverTime will be the difference in classification
%                                   between two conditions.
%       stats.StdError:             NxM matrix; standard-deviation across subjects over time-time
%       stats.pVals:                NxM matrix; p-values of each tested time-time point
%       stats.pStruct:              struct; cluster info, if mpcompcor_method was set to
%                                   'cluster_based'
%       stats.settings:             struct; the settings grabbed from the level-1 results
%       stats.condname:             string; name of format 'condition1 - condition2'.
%       stats.cfg:                  struct; the cfg of the input
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Example usage: 
%
% mask = stats1.pVals'<1; --> the pVals matrix first needs to be rotated; it consists of significant
%                             p-values within significant clusters, and 1's for all non-signficant
%                             time(-time) points; so evaluating against <1 will give a binary mask
%                             with 1's for significant, and 0's for non-significant pixels; this
%                             binary mask can be used for comparing conditions:
%
% diffstats = adam_compare_MVPA_stats(cfg,stats1,stats2,mask);
%
% part of the ADAM toolbox, by J.J.Fahrenfort, VU, 2017/2018
% 
% See also ADAM_COMPUTE_GROUP_MVPA, ADAM_MVPA_FIRSTLEVEL, ADAM_PLOT_BDM_WEIGHTS

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