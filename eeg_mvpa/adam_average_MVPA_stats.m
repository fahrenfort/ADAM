function avstats = adam_average_MVPA_stats(cfg,stats1,stats2)
% ADAM_AVERAGE_MVPA_STATS averages two stats structures that result from adam_compute_group_MVPA or
% adam_compute_group_ERP, given that the stats have the same dimensions (e.g. two time-time
% generalization matrices of equal size and same number of participants; i.e. within-subject
% repeated measures design). Also performs statistical testing against chance.
%
% Use as:
%   avstats = adam_average_MVPA_stats(cfg,stats1,stats2);
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
% part of the ADAM toolbox, by J.J.Fahrenfort, VU, 2018
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
% unpack original cfg
if isfield(stats1,'cfg')
    v2struct(stats1(1).cfg,{'fieldNames','reduce_dims','tail','cluster_pval','indiv_pval','mpcompcor_method','trainlim','testlim','reduce_dims'});
elseif isfield(stats2,'cfg')
    v2struct(stats2(1).cfg,{'fieldNames','reduce_dims','tail','cluster_pval','indiv_pval','mpcompcor_method','trainlim','testlim','reduce_dims'});
end
v2struct(cfg);
if exist('one_two_tailed','var')
    error('The cfg.one_two_tailed field has been replaced by the cfg.tail field. Please replace cfg.one_two_tailed with cfg.tail using ''both'', ''left'' or ''right''. See help for further info.');
end

% pack cfg with defaults
nameOfStruct2Update = 'cfg';
cfg = v2struct(reduce_dims,tail,cluster_pval,indiv_pval,tail,mpcompcor_method,trainlim,testlim,reduce_dims,nameOfStruct2Update);

avstats = [];
for cStats = 1:numel(stats1)
    avstats = [avstats sub_average_MVPA_stats(cfg,stats1(cStats),stats2(cStats),mask)];
end

function avstats = sub_average_MVPA_stats(cfg,stats1,stats2,mask)
% unpack cfg
v2struct(cfg);

% compute some values
ClassTotal{1} = stats1.indivClassOverTime;
ClassTotal{2} = stats2.indivClassOverTime;
nSubj = size(ClassTotal{1},1);
averageClassTotal = (ClassTotal{1}+ClassTotal{2})/2;
avstats.ClassOverTime = shiftdim(mean(averageClassTotal));
avstats.StdError = shiftdim(std(averageClassTotal)/sqrt(nSubj));
avstats.condname = ['avg of (' stats1.condname ', ' stats2.condname ')'];
settings = stats1.settings; % assuming these are the same!
settings.measuremethod = [settings.measuremethod ' average'];

% determine chance level
if any(strcmpi(settings.measuremethod,{'hr-far','dprime','hr','far','mr','cr'})) || strcmpi(stats1.cfg.plot_model,'FEM')
    chance = 0;
elseif strcmpi(settings.measuremethod,'AUC')
    chance = .5;
else
    chance = 1/settings.nconds;
end

% mask size
if isempty(mask)
	mask = ones([size(ClassTotal{1},2) size(ClassTotal{1},3)]);
end

if strcmpi(mpcompcor_method,'fdr')
    % FDR CORRECTION
    [~,ClassPvals] = ttest(averageClassTotal,chance,indiv_pval,'tail',tail);
    h = fdr_bh(squeeze(ClassPvals),cluster_pval);
    ClassPvals(~h) = 1;
    pStruct = []; % still need to implement
elseif strcmpi(mpcompcor_method,'cluster_based')
    % CLUSTER BASED CORRECTION
    [ ClassPvals, pStruct ] = cluster_based_permutation(averageClassTotal,chance,cfg,settings,mask);
elseif strcmpi(mpcompcor_method,'uncorrected')
    % NO MP CORRECTION
    [~,ClassPvals] = ttest(averageClassTotal,chance,'tail',tail);
    ClassPvals = squeeze(ClassPvals);
    ClassPvals(~mask) = 1;
    pStruct = []; % still need to implement
else
    % NO TESTING, PLOT ALL
    ClassPvals = zeros([size(ClassTotal{1},2) size(ClassTotal{1},3)]);
    pStruct = [];
end

% output stats
avstats.pVals = ClassPvals;
avstats.pStruct = pStruct;
avstats.cfg = cfg;
avstats.settings = settings;