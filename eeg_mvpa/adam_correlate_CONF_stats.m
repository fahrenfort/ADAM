function stat = adam_correlate_CONF_stats(cfg, stat1, stat2)
% adam_correlate_CONF_stats computes the correlation between the trial-by-trial, sample-by-sample
% confidence scores of the classifier stored in stat1 and either correlates these with the
% confidence scores of the classifier in another stats variable, or with some other trial by trial
% score metric (behavioral or otherwise). The resulting correlations are fisher transformed prior to
% t-testing against 0, using multiple comparison correction if desired (as in
% adam_compute_group_MVPA). The output is a stat variable that can be plotted using adam_plot_MVPA.
% 
% Usage: stat = adam_correlate_CONF_stats(cfg, stat1, stat2)
%
% Inputs:
%       cfg             A cfg structure containg the following optional fields:
%                       -   cfg.corr_method = 'Pearson' (default) to compute Pearson's linear
%                           correlation coefficient, 'Kendall' to compute Kendall's
%                           tau, or 'Spearman' to compute Spearman's rho. 
%                       -   cfg.mpcompcor_method and related fields (see adam_compute_group_MVPA)
%       stat1           A stat variable computed by adam_compute_group_MVPA using 
%                       cfg.read_confidence = true. This only works when cfg.save_confidence = true;
%                       when running adam_MVPA_firstlevel.
%       stat2           A stat variable as above, or a struct array containing the following fields:
%                       - stat2(cSubjects).scores(cTrials)      containing the to-be-correlated
%                                                               values 
%                       - stat2(cSubjects).trial_index(cTrials) containing the original trial
%                                                               indices from which these values
%                                                               originate
%                       - stat2(cSubjects).dimord               = 'rpt';
%                       - stat2(cSubjects).event_labels         contains an array with the original
%                                                               event labels of the trials. Not
%                                                               obligatory, but strongly recommended
%                                                               as a double check to ascertain that
%                                                               trials match up.
% Output:
%       stat            A group stat variable that can be plotted using adam_plot_MVPA.
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Example:
%       cfg                     = [];
%       cfg.corr_method         = 'Spearman';
%       cfg.mpcompcor_method    = 'cluster-based';
%       outstat = adam_correlate_CONF_stats(cfg,mvpa_stats(1),mvpa_stats(2));
%
% part of the ADAM toolbox, by J.J.Fahrenfort, UvA/VU, 2018
% 
% See also: ADAM_MVPA_FIRSTLEVEL, ADAM_COMPUTE_GROUP_MVPA, ADAM_PLOT_MVPA

% set defaults and unpack cfg
corr_method =       'Pearson';
mpcompcor_method =  'uncorrected';
tail =              'both';
indiv_pval =        .05;
cluster_pval =      .05;
v2struct(cfg);

% pack cfg with defaults
nameOfStruct2Update = 'cfg';
cfg = v2struct(corr_method,mpcompcor_method,tail,indiv_pval,cluster_pval,nameOfStruct2Update);

% some checks, initialize conf1 and conf2
settings = stat1.settings;
if numel(stat1)>1
    error('stat1 should be a single stat variable, not an array of stat variables.');
end
if ~isfield(stat1,'indivConf')
    error('ERROR: The first stats variable does not contain confidence scores, see the help of this function for how to resolve.');
else
    conf1 = stat1.indivConf;
end
if ~isfield(stat2,'indivConf')
    if ~isfield(stat2,'scores')
        error('ERROR: The second stats variable does not contain confidence scores, see the help of this function for how to resolve.');
    else
        conf2 = stat2;
    end
else
    conf2 = stat2.indivConf;
end
if numel(conf1) ~= numel(conf2)
    error('ERROR: the two stats variables have an unequal number of subjects.');
end
nSubj = numel(conf1);

% compute correlations
for cSubj = 1:nSubj
    indivClassOverTime(cSubj,:,:,:,:) = compute_correlations(cfg,conf1(cSubj),conf2(cSubj));
end

% perform statistics
StdError = shiftdim(squeeze(std(indivClassOverTime,0,1)/sqrt(size(indivClassOverTime,1))));
if sum(sum(StdError)) == 0 StdError = []; end % don't plot stderror when there is none
ClassOverTime = shiftdim(squeeze(mean(indivClassOverTime,1)));
chance = zeros(size(indivClassOverTime));

% fisher transform all correlations
fishTransformed = 0.5.*log((1+indivClassOverTime)./(1-indivClassOverTime));

% statistical testing
if nSubj > 1
    if strcmpi(mpcompcor_method,'fdr')
        % FDR CORRECTION
        [~,pVals] = ttest(fishTransformed,chance,indiv_pval,tail); 
        pVals = squeeze(pVals);
        h = fdr_bh(pVals,cluster_pval,'dep');
        pVals(~h) = 1;
        pStruct = compute_pstructs(h,pVals,fishTransformed,chance,cfg,settings);
    elseif strcmpi(mpcompcor_method,'cluster_based')
        % CLUSTER BASED CORRECTION
        [pVals, pStruct] = cluster_based_permutation(fishTransformed,chance,cfg,settings);
    elseif strcmpi(mpcompcor_method,'uncorrected')
        % NO MP CORRECTION
        [h,pVals] = ttest(fishTransformed,chance,indiv_pval,tail);
        pVals = squeeze(pVals);
        pStruct = compute_pstructs(squeeze(h),pVals,fishTransformed,chance,cfg,settings);
    else
        % NO TESTING, PLOT ALL
        pVals = zeros([size(indivClassOverTime,2) size(indivClassOverTime,3)]);
    end
else
    pVals = zeros([size(indivClassOverTime,2) size(indivClassOverTime,3)]);
end
pVals = shiftdim(squeeze(pVals));

% output data
stat = stat1;
stat.settings.measuremethod = [corr_method ' correlation'];
stat.settings.chance = 0;
stat.ClassOverTime = ClassOverTime;
stat.indivClassOverTime = indivClassOverTime;
stat.StdError = StdError;
stat.pVals = pVals;
stat.mpcompcor_method = mpcompcor_method;
stat.pStruct = pStruct;
stat.condname = ['corr(' stat1.condname ',' stat2.condname ')'];

function outcorr = compute_correlations(cfg,conf1,conf2)
% initialize
v2struct(cfg);

% get trial scores
scores1 = conf1.scores;
scores2 = conf2.scores;

% determine dimensions
dims1 = regexp(conf1.dimord, '_', 'split');
dims2 = regexp(conf2.dimord, '_', 'split');

% check whether trials were the first dimenions
if ~strcmpi(dims1,'rpt')
    error('Trials (rpt) should always be the first dimension when correlating classifier confidence scores.');
end

% pick classifier confidence scores for the true classes (every trial gets a confidence score for
% every class, here we only extract the confidence for the actual class of that trial)
classdim = find(strcmp(dims1,'class'));
index = cell(1,numel(dims1)); index(:) = {':'};
newindex = cell(1,numel(dims1)-1); newindex(:) = {':'};
true_labels = conf1.true_class_labels;
for cRpt = 1:numel(true_labels)
    index{1} = cRpt;
    index{classdim} = true_labels(cRpt);
    newindex{1} = cRpt;
    newscores1(newindex{:}) = scores1(index{:});
end
% removing the class dimension
dims1(classdim) = [];
scores1 = newscores1;

% do the same for scores2, if it contains separate scores for different classes
classdim = find(strcmp(dims2,'class'));
if ~isempty(classdim)
    index = cell(1,numel(dims2)); index(:) = {':'};
    newindex = cell(1,numel(dims2)-1); newindex(:) = {':'};
    true_labels = conf2.true_class_labels;
    for cRpt = 1:numel(true_labels)
        index{1} = cRpt;
        index{classdim} = true_labels(cRpt);
        newindex{1} = cRpt;
        newscores2(newindex{:}) = scores2(index{:});
    end
    % and removing the class dimension
    dims2(classdim) = [];
    scores2 = newscores2;   
end

% clear junk
clearvars new* true*;

% match trial indices in both stats (remove non matching trials)
index1 = cell(1,numel(dims1)); index1(:) = {':'};
index2 = cell(1,numel(dims2)); index2(:) = {':'};
[~, index1{1}, index2{1}] = intersect(conf1.trial_index,conf2.trial_index);
scores1 = scores1(index1{:});
scores2 = scores2(index2{:});

% display message showing the degree to which they match
num_matching = numel(index1{1});
total_trials = numel(conf1.trial_index);
perc_matching = num_matching / total_trials;
disp([ num2str(num_matching) ' out of ' num2str(total_trials) ' trials in stat variable 1 match with trials in stat variable 2']);
if perc_matching < .8
    disp(['WARNING: this is only ' num2str(perc_matching*100,3) '% of the trials. Something may be wrong, double check whether the original trial indices come from the same trials in the stats variables.']);
end

% double check whether the event codes match
if isfield(conf2,'event_labels')
    if ~all(conf1.event_labels(index1{:}) == conf2.event_labels(index2{:}))
        error('For some reason, the event codes of the trials are not matching up, re-check your pipeline.');
    end
else
    disp('WARNING: You are correlating a stats variable with behavior without using event_labels in your behavioral stats structure. This is not recommended, do make sure the trial indices match up *and* that they come from the same trials.');
end

% correlate scores across trials (first dimension), for all the dimensions in the matrix (2nd to last dimension)
size_scores1 = size(scores1);
% reshape to a two-dimensional array, so we correlate every point
scores1 = reshape(scores1,[size_scores1(1) prod(size_scores1(2:end))]);
% if they are the same size, scores2 should be reshaped too
if numel(scores1)==numel(scores2)
    scores2 = reshape(scores2,[size_scores1(1) prod(size_scores1(2:end))]);
    matchScore2 = true;
else
    % in this case, scores2 is a behavioral struct generated by the experimenter (single value for every trial)
    matchScore2 = false;
end
% now do the loop
outcorr = NaN(1,size(scores1,2));
for cInd = 1:size(scores1,2)
    if matchScore2
        outcorr(cInd) = corr(scores1(:,cInd),scores2(:,cInd),'type',corr_method);
    else
        outcorr(cInd) = corr(scores1(:,cInd),scores2,'type',corr_method);
    end
end
% reshape back to whatever dimensions it had (trials are now gone of course, as this is a correlation across trials)
outcorr = reshape(outcorr,[1 size_scores1(2:end)]);
