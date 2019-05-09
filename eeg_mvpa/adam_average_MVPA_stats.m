function avstats = adam_average_MVPA_stats(cfg,varargin)
% ADAM_AVERAGE_MVPA_STATS averages two or more stats variables that result from
% adam_compute_group_MVPA or from adam_compute_group_ERP, given that the stats have the same
% dimensions (e.g. two time-time generalization matrices of equal size and same number of
% participants; i.e. within-subject repeated measures design). Also performs statistical testing
% against chance.
%
% Use as:
%   avstats = adam_average_MVPA_stats(cfg,stats1,stats2,...);
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
%       stats1, stats2, ... =  contains two or more stats variables computed by
%                              adam_compute_group_MVPA, adam_compute_group_ERP, etc. Function
%                              assumes that stats are based on the same accuracy measure and use
%                              the same statistical test.
%
% The output diffstats structure will contain the following fields:
%
%       stats.ClassOverTime:     NxM matrix; group-average classification accuracy over N
%                                testing time points and M training time points; note that if
%                                reduce_dims is specified, M will be 1, and ClassOverTime will be
%                                squeezed to a Nx1 matrix of classification over time. Here,
%                                ClassOverTime will be the average classification performance
%                                over the averaged conditions.
%       stats.indivClassOverTime SxNxM; same as above over S subjects.
%       stats.StdError:          NxM matrix; standard-error across subjects over time-time
%       stats.pVals:             NxM matrix; p-values of each tested time-time point
%       stats.pStruct:           struct; cluster info, currently only appears if mpcompcor_method
%                                was set to 'cluster_based'
%       stats.settings:          struct; the settings used during the level-1 (single subject)
%                                results
%       stats.condname:          string; name of format 'average(condition1,condition2)'.
%       stats.cfg:               struct; the cfg used to create these stats
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% part of the ADAM toolbox, by J.J.Fahrenfort, VU, 2018
%
% See also ADAM_COMPUTE_GROUP_MVPA, ADAM_MVPA_FIRSTLEVEL, ADAM_PLOT_MVPA, ADAM_PLOT_BDM_WEIGHTS, FDR_BH

if nargin<2
    disp('cannot plot graph without stats input, need at least 2 arguments:');
    help adam_plot_MVPA;
    return
end

% concatenate stats
stats = concat_stats(varargin{:});

% get some defaults
plot_model = 'BDM';
mask = [];
reduce_dims = '';
tail = 'both';
cluster_pval = .05;
indiv_pval = .05;
mpcompcor_method = 'uncorrected';
plot_dim = 'time_time';
trainlim = [];
testlim = [];
% unpack original cfg
if isfield(stats(1),'cfg')
    v2struct(stats(1).cfg);
end
v2struct(cfg);
if exist('one_two_tailed','var')
    error('The cfg.one_two_tailed field has been replaced by the cfg.tail field. Please replace cfg.one_two_tailed with cfg.tail using ''both'', ''left'' or ''right''. See help for further info.');
end

% pack cfg with defaults
nameOfStruct2Update = 'cfg';
cfg = v2struct(reduce_dims,tail,cluster_pval,indiv_pval,tail,mpcompcor_method,trainlim,testlim,reduce_dims,plot_dim,nameOfStruct2Update);

% unpack cfg
v2struct(cfg);
nStats = numel(stats);
% initialize decoding vars
indivClassOverTime = zeros(size(stats(1).indivClassOverTime));
ClassOverTime = zeros(size(stats(1).ClassOverTime));
% initialize weight vars
getweights = isfield(stats(1),'weights') && stats(1).settings.FEM;
if getweights
    avWeights = zeros(size(stats(1).weights.avWeights));
    indivWeights = zeros(size(stats(1).weights.indivWeights));
    CTF = zeros(size(stats(1).weights.CTF));
    indivCTF = zeros(size(stats(1).weights.indivCTF));
    nConds = numel(stats(1).weights.CTFpercond);
    for cCond = 1:nConds
        CTFpercond{cCond} = zeros(size(stats(1).weights.CTFpercond{cCond}));
        indivCTFpercond{cCond} = zeros(size(stats(1).weights.indivCTFpercond{cCond}));
    end
end
% add up to compute average of stats
for cStats = 1:nStats
    ClassOverTime = ClassOverTime + stats(cStats).ClassOverTime;
    indivClassOverTime = indivClassOverTime + stats(cStats).indivClassOverTime;
    % also compute average of weights
    if getweights
        avWeights = avWeights + stats(cStats).weights.avWeights;
        indivWeights = indivWeights + stats(cStats).weights.indivWeights;
        CTF = CTF + stats(cStats).weights.CTF;
        indivCTF = indivCTF + stats(cStats).weights.indivCTF;
        for cCond = 1:nConds
            CTFpercond{cCond} = CTFpercond{cCond} + stats(cStats).weights.CTFpercond{cCond};
            indivCTFpercond{cCond} = indivCTFpercond{cCond} + stats(cStats).weights.indivCTFpercond{cCond};
        end
    end
end
% average and store
nSubj = size(indivClassOverTime,1);
ClassOverTime = ClassOverTime/nStats;
indivClassOverTime = indivClassOverTime/nStats;
avstats.ClassOverTime = ClassOverTime;
avstats.indivClassOverTime = indivClassOverTime;
if getweights
    avstats.weights.avWeights = avWeights/nStats;
    avstats.weights.indivWeights = indivWeights/nStats;
    avstats.weights.CTF = CTF/nStats;
    avstats.weights.indivCTF = indivCTF/nStats;
    for cCond = 1:nConds;
        avstats.weights.CTFpercond{cCond} = CTFpercond{cCond}/nStats;
        avstats.weights.indivCTFpercond{cCond} = indivCTFpercond{cCond}/nStats;
    end
end
% re-compute SEMs
if nSubj > 1
    avstats.StdError = shiftdim(std(avstats.indivClassOverTime)/sqrt(nSubj));
    if getweights
        avstats.weights.semCTF = shiftdim(std(avstats.weights.indivCTF))/sqrt(nSubj);
        for cCond = 1:nConds
           avstats.weights.semCTFpercond{cCond} = shiftdim(std(avstats.weights.indivCTFpercond{cCond}))/sqrt(nSubj);
        end
    end
else
    avstats.StdError = [];
end

avstats.condname = ['average(' cell2csv({stats(:).condname}) ')'];
settings = stats(1).settings; % assuming these are the same!

% determine chance level
if ~isfield(settings,'chance')
    if any(strcmpi(settings.measuremethod,{'hr-far','dprime','hr','far','mr','cr'})) || strncmpi(settings.measuremethod,'\muV',4) || ~isempty(strfind(settings.measuremethod,'difference')) || strcmpi(plot_model,'FEM')
        chance = 0;
    elseif strcmpi(settings.measuremethod,'AUC')
        chance = .5;
    else
        chance = 1/settings.nconds;
    end
    settings.chance = chance;
else
    chance = settings.chance;
end

% mask size
if isempty(mask)
    mask = ones([size(indivClassOverTime,2) size(indivClassOverTime,3)]);
end

if nSubj > 1
    if strcmpi(mpcompcor_method,'fdr')
        % FDR CORRECTION
        [~,ClassPvals] = ttest(indivClassOverTime,chance,indiv_pval,tail);
        ClassPvals = squeeze(ClassPvals);
        h = fdr_bh(ClassPvals,cluster_pval,'dep');
        ClassPvals(~h) = 1;
        pStruct = compute_pstructs(h,ClassPvals,indivClassOverTime,chance,cfg,settings);
    elseif strcmpi(mpcompcor_method,'cluster_based')
        % CLUSTER BASED CORRECTION
        [ ClassPvals, pStruct ] = cluster_based_permutation(indivClassOverTime,chance,cfg,settings,mask);
    elseif strcmpi(mpcompcor_method,'uncorrected')
        % NO MP CORRECTION
        [h, ClassPvals] = ttest(indivClassOverTime,chance,indiv_pval,tail);
        ClassPvals = squeeze(ClassPvals);
        pStruct = compute_pstructs(squeeze(h),ClassPvals,indivClassOverTime,chance,cfg,settings);
    else
        % NO TESTING, PLOT ALL
        ClassPvals = zeros([size(indivClassOverTime,2) size(indivClassOverTime,3)]);
    end
else
    ClassPvals = zeros([size(indivClassOverTime,2) size(indivClassOverTime,3)]);
end
ClassPvals = shiftdim(squeeze(ClassPvals));

% output stats
avstats.pVals = ClassPvals;
if exist('pStruct','var')
    avstats.pStruct = pStruct;
end
% extract chancel level in order to be able to compute latencies (i.e. subtract chance during latency computation)
if ~isfield(settings,'chance') % backwards compatibility
    if any(strcmpi(settings.measuremethod,{'hr-far','dprime','hr','far','mr','cr'})) || strcmpi(settings.measuremethod,'\muV') || ~isempty(strfind(settings.measuremethod,' difference')) || ~isempty(strfind(settings.measuremethod,' correlation')) || strcmpi(plot_model,'FEM')
        cfg.chance = 0;
    elseif strcmpi(settings.measuremethod,'AUC')
        cfg.chance = .5;
    else
        cfg.chance = 1/settings.nconds;
    end
else
    cfg.chance = settings.chance;
end
avstats.cfg = cfg;
avstats.settings = settings;
% compute latency
try
    avstats.latencies = extract_latency(cfg,avstats);
catch ME
    disp('Cannot extract latencies.');
    disp(ME.message);
    avstats.latencies = [];
end