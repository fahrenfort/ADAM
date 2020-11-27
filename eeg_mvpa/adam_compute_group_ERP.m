function stats = adam_compute_group_ERP(cfg,folder_name)
% ADAM_COMPUTE_GROUP_ERP computes a group-level event-related potential (ERP) from single-subject
% ERPs that were computed by adam_MVPA_firstlevel. It also performs statistical tests at the group
% level. When executing adam_compute_group_ERP, a selection dialog will pop up. This dialog allows
% the user to select a directory containing the single subject results (output of
% adam_MVPA_firstlevel) for which to compute the group stats variable. The user can either select a
% directory referring to a specific analysis (e.g. EEG_FAM_VS_NONFAMOUS), or select one directory
% higher up in the hierarchy (e.g. RAW_EEG) which may contains several such analyses. The function
% creates a stats output variable. This variable is a structure that contains the group average,
% together with the single subject ERPs, and (if requested) statistical outcomes, and can be used as
% input for the adam_plot_MVPA plotting function. When selecting a directory containing several
% analyses, the stats variable will be an array of which each element contains a group ERP.
% adam_compute_group_ERP has several options to extract ERPs belonging to particular conditions or
% from particular electrodes, average or subtract signals from different electrodes etc. These
% options are outlined below.
%
% Use as:
%   stats = adam_compute_group_ERP(cfg)
%
% The cfg (configuration) input structure should specify electrode selection and statistics options.
% The cfg can contain the following optional parameters:
%
%
%       cfg.startdir         = '' (default); ADAM will pop-up a selection dialog when running
%                              adam_compute_group_MVPA. The cfg.startdir parameter allows you to
%                              specify the starting directory of this selection dialog. Use this
%                              parameter to specify where the results of all the first level
%                              analyses are located. When you do not specify cfg.startdir, you will
%                              be required to navigate from your Matlab root folder to the desired
%                              results directory every time you run a group analysis.
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
%       cfg.mask            =  Optionally, you can provide a mask: a vector to pre-select a 'region
%                              of interest' to constrain the comparison. You can for example base
%                              the mask on pVals matrix in stats.pVals (beware of double dipping).
%       cfg.electrode_def    = cell array of strings specifying which electrodes to extract, e.g.
%                              {'O1','Oz','O2'}, or a cell array of cell arrays, e.g.:
%                              {{'PO7'},{'PO8'};{'PO8'},{'PO7'}} for a lateralization analysis
%                              (examples below); specifying one electrode or a group of
%                              electrodes that are then averaged or subracted; make sure the
%                              labeling you specify corresponds to the labels present in the raw
%                              data that you used as input for adam_MVPA_firstlevel.
%       cfg.electrode_method = 'average' (default); string specifying what to do if multiple
%                              electrodes are specificied: 'average', the ERP of the average of the
%                              specified electrodes is computed; 'subtract', two sets of electrodes
%                              or two invidual electrodes are first averaged and then subtracted,
%                              e.g. if you want to do a lateralization analysis (left versus right
%                              parietal-occipital channels); with one electrode, there is obviously
%                              nothing to average or subtract, so 'average' is the default).
%       cfg.condition_def    = vector of integers specifying which conditions (stimulus classes
%                              that served as input to the first-level analysis) to extract,
%                              substract (compare), or or average; e.g. [1,2,3,4].
%       cfg.condition_method = 'keep' (default); string specifying what to do with the requested
%                              conditions; in case of 'keep', the ERP of each condition (stimulus
%                              class of the first-level analysis) is computed separately and tested
%                              against zero; other options are 'average' wich performs a
%                              condition-average against zero, or 'subtract' in case of two
%                              conditions specified in condition_def, which performs a
%                              condition-comparison: the second condition is first subtracted from
%                              the first condition, and the result is tested against zero.
%       cfg.timelim          = [min max]; vector specifiyin a limited time range to analyze, if
%                              desired; if not specified, the whole time range is used.
%       cfg.plotsubjects     = false (default); or true; if true, the individual ERPs of all
%                              subjects that are extracted will be plotted in a single figure, with
%                              a separate subplot for each subject to enable inspection of the data
%                              underlying group averages.
%       cfg.resample_eeg     = integer for a new down-sampled sampling rate if desired; default: 0
%                              (no resampling).
%
% The output stats structure will contain the following fields:
%
%       stats.ClassOverTime:    1xN matrix; group-average ERP over N time points (name ClassOverTime
%                               comes from MVPA nomenclature: Classification over time).
%       stats.StdError:         1xN matrix; standard-deviation across subjects over time
%       stats.pVals:            1xN matrix; p-values of each tested time point
%       stats.pStruct:          struct; cluster info, if mpcompcor_method was set to
%                               'cluster_based'
%       stats.latencies         struct; containing three fields with onset latencies, based on
%                               the percent-amplitude latency (50%): single subject, jackknife and
%                               GA (grand average). Note that single subject outcomes are not very
%                               reliable, it is much better to rely on the jackknife results. For
%                               details about latency computation see: Liesefeld, H. R. (2018).
%                               Estimating the Timing of Cognitive Operations With MEG/EEG Latency
%                               Measures: A Primer, a Brief Tutorial, and an Implementation of
%                               Various Methods. Frontiers in Neuroscience, 12, 765.
%                               http://doi.org/10.3389/fnins.2018.00765
%       stats.mpcompcor_method: string; correction method ('uncorrected' is default)
%       stats.settings:         struct; the settings grabbed from the level-1 results
%       stats.condname:         string; combining name of the level-1 folder and the condition_method
%       stats.channelpool:      string; summarizing the specified electrodes
%       stats.reduce_dims:      [] (only relevant for MVPA)
%       stats.cfg:              struct; the cfg used to create these stats
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Common use cases:
%   (1) get ERP difference between conditions
%   (2) get condition-average ERP and test against zero
%   (3) get condition-specific ERPs and test each against zero
%   (4) get N2pc/CDA electrode subtractions for two conditions, average them, and test against zero
%
% Use case (1):
% cfg.electrode_def    = { 'Oz', 'Iz', 'POz' };
% cfg.elecrode_method  = 'average' (default)
% cfg.condition_def    = [1,2];
% cfg.condition_method = 'subtract'
%
% Use case (2):
% cfg.electrode_def    = { 'Oz' };
% cfg.condition_def    = [1,2,3,4];
% cfg.condition_method = 'average'
%
% Use case (3):
% cfg.electrode_def    = { 'Oz' };
% cfg.condition_def    = [1,2,3,4];
% cfg.condition_method = 'keep' (default)
%
% Use case (4):
% cfg.condition_def     = [1,2];
% cfg.condition_methods = 'average'
% cfg.electrode_def     = {{'PO7'},{'PO8'};{'PO8'},{'PO7'}};
% cfg.electrode_method  = 'subtract'; -->   specifies subtraction method, such that for condition 1
%                                           PO8 is subtracted from P07, while for condition 2 P07 is
%                                           subtracted from PO8. These are subsequently averaged.
%
% part of the ADAM toolbox, by J.J.Fahrenfort, VU, 2017/2018
%
% See also ADAM_MVPA_FIRSTLEVEL, ADAM_COMPUTE_GROUP_MVPA, ADAM_PLOT_MVPA, ADAM_PLOT_BDM_WEIGHTS, FDR_BH

if nargin<2
    folder_name = '';
end
plot_order = {};

% backwards compatibility
v2struct(cfg);
if exist('one_two_tailed','var')
    error('The cfg.one_two_tailed field has been replaced by the cfg.tail field. Please replace cfg.one_two_tailed with cfg.tail using ''both'', ''left'' or ''right''. See help for further info.');
end
if exist('plotorder','var')
    plot_order = plotorder;
    cfg.plot_order = plot_order;
    cfg = rmfield(cfg,'plotorder');
end
if exist('channels','var') && ~isfield(cfg,'channelpool')
    channelpool = channels;
    cfg.channelpool = channelpool;
end

% Main routine, is a folder name specified? If not, pop up selection dialog
if isempty(folder_name)
    if ~isfield(cfg,'startdir')
        cfg.startdir = '';
        disp('NOTE: it is easier to select a directory when you indicate a starting directory using cfg.startdir, otherwise you have to start selection from root every time...');
    end
    folder_name = uigetdir(cfg.startdir,'select directory for which to compute group ERP results');
    if ~ischar(folder_name)
        error('no folder was selected');
    end
end
if ~exist(folder_name,'dir')
    error('the specified folder does not exist');
end
cfg.folder = folder_name;

% where am I?
ndirs = drill2data(folder_name);
if isempty(plot_order)
    dirz = dir(folder_name);
    dirz = {dirz([dirz(:).isdir]).name};
    plot_order = dirz(~(strcmp(dirz,'.')|strcmp(dirz,'..')));
    if ndirs == 1
        [folder_name, plot_order] = fileparts(folder_name);
        plot_order = {plot_order};
    elseif ndirs > 2
        error('You seem to be selecting a directory that is too high in the hiearchy, drill down a little more.');
    end
    cfg.plot_order = plot_order;
elseif ndirs ~= 2
    error('You seem to be selecting a directory that is either too high or too low in the hiearchy given that you have specified cfg.plot_order. Either remove cfg.plot_order or select the appropriate level in the hierarchy.');
else
    dirz = dir(folder_name);
    dirz = {dirz([dirz(:).isdir]).name};
    dirz = dirz(~(strcmp(dirz,'.')|strcmp(dirz,'..')));
    for cPlot = 1:numel(plot_order)
        dirindex = find(strcmpi(plot_order{cPlot},dirz));
        if isempty(dirindex) % if an exact match cannot be made, look only for the pattern in the first sequency of characters
            dirindex = find(strcmpi(plot_order{cPlot},dirz,numel(plot_order{cPlot})));
        end
        if isempty(dirindex)
            error(['cannot find condition ' plot_order{cPlot} ' specified in cfg.plot_order']);
        elseif numel(dirindex) > 1
            error(['cannot find a unique condition for the pattern ' plot_order{cPlot} ' specified in cfg.plot_order']);
        else
            plot_order{cPlot} = dirz{dirindex};
        end
    end
    if ~all(ismember(plot_order,dirz))
        error('One or more of the folders specified in cfg.plot_order cannot be found in this results directory. Change cfg.plot_order or select a different directory for plotting.');
    end
end

% loop through directories (results folders)
stats = [];
for cdirz = 1:numel(plot_order)
    stats = [stats subcompute_group_ERP(cfg,[folder_name filesep plot_order{cdirz}])];
end

% fill cfg.plot_order in case not given by user and from single folder
if numel(cfg.plot_order) == 1
    cfg.plot_order = {stats(:).condname};
end

function [stats,cfg] = subcompute_group_ERP(cfg,folder_name)
% set defaults
tail = 'both';
indiv_pval = .05;
cluster_pval = .05;
plotsubjects = false;
name = [];
reduce_dims = [];
mpcompcor_method = 'uncorrected';
timelim = [];
resample_eeg = 0;
electrode_def = [];
condition_def = [1,2]; % By default substracting cond1 - cond2
electrode_method = 'keep';
condition_method = 'keep';
% unpack graphsettings
plottype = '2D';
channelpool = '';
exclsubj = [];
v2struct(cfg);

% pack graphsettings with defaults
nameOfStruct2Update = 'cfg';
cfg = v2struct(tail,reduce_dims,indiv_pval,cluster_pval,name,plottype,mpcompcor_method,electrode_def,electrode_method,condition_def,condition_method,timelim,resample_eeg,nameOfStruct2Update);

% fill some empties
if isempty(electrode_def)
    error('no electrode_def was specified in cfg, set cfg.electrode_def to some electrode, e.g. ''Oz'', see help of this function for more info.');
end
if isempty(channelpool)
    chandirz = dir(folder_name);
    chandirz = {chandirz([chandirz(:).isdir]).name};
    chandirz = sort(chandirz(cellfun(@isempty,strfind(chandirz,'.'))));
    channelpool = chandirz{1};
    cfg.channelpool = channelpool;
    disp(['No cfg.channelpool specified, defaulting to channelpool ' channelpool ]);
end

% get filenames
plotFreq = ''; % this is empty for now, but might be used to look at the ERPs obtained from a TF analysis
subjectfiles = dir([folder_name filesep channelpool plotFreq filesep '*.mat']);
[~, condname,ext] = fileparts(folder_name);
condname = [condname ext];
subjectfiles = { subjectfiles(:).name };
subjectfiles(strncmp(subjectfiles,'.',1)) = []; % remove hidden files

% limiting subjects
if ~isempty(exclsubj)
    subjectfiles = select_subjects(subjectfiles,exclsubj,true);
end

% see if data exists
nSubj = numel(subjectfiles);
if nSubj == 0
    error(['cannot find data in specified folder ' folder_name filesep channelpool plotFreq ' maybe you should specify (a different) cfg.channelpool?']);
end

% prepare figure in case individual subjects are plotted
if plotsubjects
    fh = figure('name',['individual subjects, condition: ' condname]);
    set(fh, 'Position', get(0,'Screensize'));
    set(fh,'color','w');
end

% do the loop, restrict time and frequency if applicable
for cSubj = 1:nSubj
    fprintf(1,'loading subject %d of %d\n', cSubj, nSubj);
    matObj = matfile([folder_name filesep channelpool plotFreq filesep subjectfiles{cSubj}]);
    if ~isempty(whos(matObj,'settings'))
        settings = matObj.settings;
    else % BW compatibility, will become obsolete over time
        error('Cannot find settings struct');
    end
    % OLD for backwards compatible
    if isfield(settings,'ft_erpstruct')
        FT_ERP = settings.ft_erpstruct;
    elseif isfield(settings,'FT_ERP')
        FT_ERP = settings.FT_ERP;
    elseif ~isempty(whos(matObj,'FT_ERP'))
        FT_ERP = matObj.FT_ERP;
    else
        error('no ERP was computed for these data');
    end
    
    if iscell(FT_ERP)
        FT_ERP = FT_ERP{2};
        if isfield(FT_ERP,'origindex')
            FT_ERP = rmfield(FT_ERP,'origindex');
        end
        if isfield(FT_ERP,'oldindex')
            FT_ERP = rmfield(FT_ERP,'oldindex');
        end
        disp('Extracting ERPs from testing data');
    end
    
    % extract electrodes and/or do electrode subtraction/averaging
    FT_ERP = restrict_FT_ERP(FT_ERP,cfg);
    if numel(FT_ERP) == 1
        settings.times = {FT_ERP.time, FT_ERP.time};
    else
        settings.times = {FT_ERP{1}.time, FT_ERP{1}.time};
    end
    
    % possible use cases:
    % 	(1) get ERP difference between conditions + stats against zero for one or more results folders
    % 	(2) get raw ERPs for one or more conditions from one or more results folders and average them
    %   (3) get raw ERPs for more conditions from a given results folder without averaging
    % 	(4) get N2pc/CDA electrode subtractions for two conditions, average
    %       them, and do stats against zero (for each results folder)
    %   (5) extract raw ERPs for multiple electrodes of a single condition
    
    % now do condition subtraction
    clear trial;
    if strcmpi(condition_method,'subtract')
        % only works for a single electrode
        if numel(FT_ERP)>1
            error('can only subtract conditions from single electrode, not multiple electrodes');
        end
        % subtract
        if ~(size(condition_def,2)==2)
            error('Condition_def does not contain the correct number of conditions (2) to be able to subtract.');
        end
        trial = FT_ERP.trial(FT_ERP.trialinfo==condition_def(1),:,:) - FT_ERP.trial(FT_ERP.trialinfo==condition_def(2),:,:);
    elseif strcmpi(condition_method,'average')
        % only works for a single electrode
        if numel(FT_ERP)>1
            error('can only average conditions from single electrode, not multiple electrodes');
        end
        % average
        trial = mean(FT_ERP.trial(ismember(FT_ERP.trialinfo,condition_def),:,:),1);
    else
        if numel(FT_ERP) == 1
            % plain extraction of all conditions
            trial = FT_ERP.trial(ismember(FT_ERP.trialinfo,condition_def),:,:);
        else
            % let's pretend that electrodes are conditions so we can get multiple electrodes out
            for cEl = 1:numel(FT_ERP)
                trial(cEl,:) = FT_ERP{cEl}.trial(ismember(FT_ERP{cEl}.trialinfo,condition_def),:,:);
            end
        end
    end
    
    % little hack to change V to muV (if applicable)
    if  all(all(all(trial<10^-3)))
        trial = trial*10^6;
    end
    
    % hold on to each subject's resulting erps! A bit hackey, but the cCond can also refer to electrodes, 
    % if multiple electrodes for a single condition were extracted
    for cCond=1:size(trial,1)
        ClassTotal{cCond}(cSubj,:) = squeeze(trial(cCond,:));
    end
    
    % plot individual subjects
    if plotsubjects
        subplot(numSubplots(nSubj,1),numSubplots(nSubj,2),cSubj);
        onestat.ClassOverTime = squeeze(trial);
        onestat.StdError = [];
        onestat.pVals = ones(size(squeeze(trial)));
        onestat.indivClassOverTime = [];
        onestat.settings = settings;
        onestat.settings.measuremethod = '\muV';
        onestat.condname = condname;
        onestat.channelpool = channelpool;
        onestat.cfg = [];
        tmpcfg = cfg;
        if isempty(matObj.BDM)
            plot_model = 'FEM';
        else
            plot_model = 'BDM';
        end
        tmpcfg = cfg;
        tmpcfg.plot_model = plot_model;
        tmpcfg.plotsubjects = true;
        tmpcfg.plot_order = {condname};
        tmpcfg.acclim2D = [];
        tmpcfg.acclim3D = [];
        tmpcfg.acctick = [];
        tmpcfg.ploterp = true;
        adam_plot_MVPA(tmpcfg,onestat);
        subjname = subjectfiles{cSubj};
        underscores = strfind(subjname,'_');
        subjname = regexprep(subjname(underscores(2)+1:underscores(end)-1),'_',' ');
        ntitle(subjname,'fontsize',10,'fontweight','bold');
        drawnow;
    end
end

% next, compute stats
chance = 0;

% statistical testing
origcondname = condname;
for cCond = 1:numel(ClassTotal) % loop over stats

    % determine condname
    if strcmpi(electrode_method,'subtract') && strcmpi(condition_method,'average')
        condname =  [origcondname '-average(' cell2csv(FT_ERP.channelpool) ')'];
    elseif strcmpi(electrode_method,'subtract')
        condname =  [origcondname ' erp' num2str(condition_def(cCond)) '-' cell2csv(FT_ERP.channelpool)];
    elseif strcmpi(condition_method,'subtract')
        condname =  [origcondname '-subtraction(' cell2csv(FT_ERP.channelpool) ')'];
    elseif strcmpi(condition_method,'average')
        condname = [origcondname ' average(' cell2csv(FT_ERP.channelpool) ')'];
    else
        if numel(FT_ERP) == 1
            condname = [origcondname ' erp' num2str(condition_def(cCond)) ' (' cell2csv(FT_ERP.channelpool) ')'];
        else % little hack for when multiple electrodes of a single condition were extracted:
            condname = [origcondname ' erp' num2str(condition_def(1)) ' (' cell2csv(FT_ERP{cCond}.channelpool) ')'];
        end
    end
    
    % get some stats
    indivClassOverTime = ClassTotal{cCond};
    ClassAverage = shiftdim(squeeze(mean(ClassTotal{cCond},1)));
    ClassStdErr = shiftdim(squeeze(std(ClassTotal{cCond},0,1)/sqrt(size(ClassTotal{cCond},1))));
    
    if strcmp(mpcompcor_method,'fdr')
        % FDR CORRECTION
        [~,ClassPvals] = ttest(ClassTotal{cCond},chance,indiv_pval,tail);
        ClassPvals = squeeze(ClassPvals);
        h = fdr_bh(ClassPvals,cluster_pval,'dep');
        ClassPvals(~h) = 1;
        pStruct = compute_pstructs(h,ClassPvals,ClassTotal{cCond},chance,cfg,settings);
    elseif strcmp(mpcompcor_method,'cluster_based')
        % CLUSTER BASED CORRECTION
        [ClassPvals, pStruct] = cluster_based_permutation(ClassTotal{cCond},chance,cfg,settings);
        % compute Pstruct
    elseif strcmp(mpcompcor_method,'uncorrected')
        % NO MP CORRECTION
        [h,ClassPvals] = ttest(ClassTotal{cCond},chance,indiv_pval,tail);
        ClassPvals = squeeze(ClassPvals);
        pStruct = compute_pstructs(h,ClassPvals,ClassTotal{cCond},chance,cfg,settings);
    else
        % NO TESTING, PLOT ALL
        ClassPvals = zeros(1,size(ClassTotal{cCond},2));
        pStruct = [];
    end
    ClassPvals = shiftdim(squeeze(ClassPvals));
    
    % outputs: put it in a matrix for consistency in plot function
    settings.measuremethod = '\muV';
    settings.chance = 0;
    stats(cCond).ClassOverTime = ClassAverage;
    stats(cCond).indivClassOverTime = indivClassOverTime;
    stats(cCond).StdError = ClassStdErr;
    stats(cCond).pVals = ClassPvals;
    stats(cCond).mpcompcor_method = mpcompcor_method;
    stats(cCond).settings = settings;
    stats(cCond).condname = condname;
    if numel(FT_ERP) == 1
        stats(cCond).channelpool = FT_ERP.channelpool;
    else % little hack for when multiple electrodes of a single condition were extracted:
        stats(cCond).channelpool = FT_ERP{cCond}.channelpool;
    end
    stats(cCond).pStruct = pStruct;
    % compute latency
    try
        stats(cCond).latencies = extract_latency(cfg,stats(cCond));
    catch ME
        disp('Cannot extract latencies.');
        disp(ME.message);
        stats(cCond).latencies = [];
    end
    stats(cCond).reduce_dims = reduce_dims;
    stats(cCond).cfg = cfg;
    if isfield(stats(cCond).cfg,'plotsubjects')
        stats(cCond).cfg = rmfield(stats(cCond).cfg,'plotsubjects');
    end
end


disp('done!');

function [FT_EEG] = restrict_FT_ERP(FT_EEG,cfg)
% resample / restrict the ERP
electrode_method = 'average';
v2struct(cfg);
% resample?
if resample_eeg
    cfg = [];
    cfg.resamplefs = resample_eeg;
    FT_EEG = ft_resampledata(cfg,FT_EEG);
end
% limit time?
if ~isempty(timelim)
    FT_EEG = select_time_from_FT_EEG(FT_EEG,(FT_EEG.time>min(timelim)/1000 & FT_EEG.time<max(timelim)/1000)); % time should be in seconds in FT_EEG
end
clear channelpool;
if strcmpi(electrode_method,'subtract')
    % subtracting electrode sets, subtracts electrode 2 from electrode 1 for each condition
    if size(electrode_def,2)~=2
        error('to subtract electrode sets, you should define two columns in cfg.electrode_def');
    end
    if size(electrode_def,1) ~= numel(condition_def)
        electrode_def = repmat(electrode_def,[numel(condition_def),1]);
    end
    for cCond = 1:size(electrode_def,1)
        % first extracting specified electrodes for specific conditions
        for cDif=1:size(electrode_def,2)
            cfg = [];
            FT_TEMP(cDif) = select_channels_from_FT_EEG(FT_EEG,electrode_def{cCond,cDif}); % this is a bare extraction
        end
        % next subtracting averages of the two electrodes for the relevant conditions
        % This was a bug: trial(cCond,:,:) = FT_TEMP(1).trial(FT_TEMP(1).trialinfo==condition_def(cCond),:,:) - FT_TEMP(2).trial(FT_TEMP(2).trialinfo==condition_def(cCond),:,:);
        % This was a bug: channelpool{cCond} = [ 'chan_subtract(' FT_TEMP(1).label ',' FT_TEMP(2).label ')'];
        trial(FT_TEMP(1).trialinfo==condition_def(cCond),:,:) = FT_TEMP(1).trial(FT_TEMP(1).trialinfo==condition_def(cCond),:,:) - FT_TEMP(2).trial(FT_TEMP(2).trialinfo==condition_def(cCond),:,:);
        channelpool{FT_TEMP(1).trialinfo==condition_def(cCond)} = [ 'chan_subtract(' FT_TEMP(1).label ',' FT_TEMP(2).label ')'];
    end
    FT_EEG.trial = trial;
    FT_EEG.channelpool = channelpool;
elseif strcmpi(electrode_method,'average')
    % averaging electrodes
    if ~iscell(electrode_def)
        electrode_def = {electrode_def};
    end
    FT_EEG = select_channels_from_FT_EEG(FT_EEG,electrode_def);
    if numel(electrode_def) > 1
        FT_EEG.channelpool = ['average(' cell2csv(FT_EEG.label) ')'];
    else
        FT_EEG.channelpool = FT_EEG.label;
    end
else
    % extracting individual electrodes (plain extract) 
    % -> only works for single condition
    if numel(condition_def) > 1
        error('You are attempting to extract multiple electrodes for multiple conditions. This is not supported. Either specify cfg.electrode_method as ''average'' (default) or ''subtract'' and/or specify cfg.condition_def for a single condition.');
    end
    if ~iscell(electrode_def)
        electrode_def = {electrode_def};
    end
    clear FT_TEMP;
    for cEl = 1:numel(electrode_def)
        % extract separate electrodes
        FT_TEMP{cEl} = select_channels_from_FT_EEG(FT_EEG,electrode_def{cEl});
        FT_TEMP{cEl}.channelpool = FT_TEMP{cEl}.label;
    end
    % pass back only the relevant electrodes as an array
    FT_EEG = FT_TEMP;
end
if numel(FT_EEG) == 1 && any(size(FT_EEG.trial)==0)
    dims = regexp(FT_EEG.dimord,'_','split');
    error(['There are no dimensions left in these fields, change cfg selection parameters for ''' cell2csv(dims(logical(size(FT_EEG.trial)==0))) '''.']);
end

function ndirs = drill2data(folder_name)
% drills down until it finds data, returns the number of directories it had
% to drill
notfound = true;
ndirs = 0;
while notfound
    dirz = dir(folder_name);
    dirz = {dirz([dirz(:).isdir]).name};
    nextlevel = dirz(cellfun(@isempty,strfind(dirz,'.')));
    if isempty(nextlevel)
        error('cannot find data, check path settings');
    end
    folder_name = fullfile(folder_name,nextlevel{1});
    ndirs = ndirs + 1;
    containsmat = ~isempty(dir(fullfile(folder_name, '*.mat')));
    containsfreq = ~isempty(dir(fullfile(folder_name, 'freq*'))) || ~isempty(dir(fullfile(folder_name, 'allfreqs')));
    if containsmat || containsfreq
        notfound = false;
    end
end

