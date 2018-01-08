function [stats,cfg] = adam_compute_group_ERP(cfg,folder_name)
% ADAM_COMPUTE_GROUP_ERP computes a group-level event-related potential (ERP) from single-subject
% ERP results from ADAM_MVPA_FIRSTLEVEL, and if requested performs a statistical test at the group
% level. The result is a 1xN structure for N conditions that contains the group average, together
% with the subject-specific ERPs and (if requested) statistical outcomes, and can be used for
% plotting (as an input for ADAM_PLOT_MVPA).
%
% Use as:
%   stats = adam_compute_group_ERP(cfg)
%
% The cfg (configuration) input structure should specify electrode selection and statistics options,
% as well as the path to the directory where the results are located. The following options can be
% specified in the cfg:
%
%
%       cfg.startdir         = string specifiying a directory where the results of 
%                              ADAM_MVPA_FIRSTLEVEL are located;
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
%       cfg.electrode_def    = string between curly brackets, e.g. {'O1','Oz','O2'}, or embedded 
%                              curly brackets: {{'PO7'},{'PO8'};{'PO8'},{'PO7'}} for a
%                              lateralization analysis; specifying one electrode or a group of
%                              electrodes that are then averaged or subracted; make sure the
%                              labeling you specify corresponds to the labels present in the raw
%                              data that you used as input for ADAM_MVPA_FIRSTLEVEL.
%       cfg.electrode_method = 'average' (default); string specifying what to do if multiple 
%                              electrodes are specificied: 'average', the ERP of the average of the
%                              specified electrodes is computed; 'subtract', two sets of electrodes
%                              or two invidual electrodes are first averaged and then subtracted,
%                              e.g. if you want to do a lateralization analysis (left versus right
%                              parietal-occipital channels); with one electrode, there is obviously
%                              nothing to average or subtract, so 'average' is the default).
%       cfg.timelim          = [min max]; vector specifiyin a limited time range to analyze, if
%                              desired; if not specified, the whole time range is used.
%       cfg.condition_def    = vector of integers specifying which conditions to take and to compare
%                              or average; e.g. [1,2,3,4].
%       cfg.condition_method = 'keep' (default); string specifying what to do with the requested
%                              conditions; in case of 'keep', the ERP of each condition is computed
%                              separately and tested against baseline; other options are 'average'
%                              wich performs a condition-average against baseline, or 'subtract' in
%                              case of two conditions specified in condition_def, which performs a
%                              condition-comparison: the second condition is first subtracted from
%                              the first condition, and the result is tested against zero.
%       cfg.plotsubjects     = false (default); or true; if true, during importing single-subject
%                              data, one figure with subplots is generated that displays each
%                              single-subject ERP.
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
%       stats.mpcompcor_method: string; correction method ('uncorrected' is default)
%       stats.settings:         struct; the settings grabbed from the level-1 results
%       stats.condname:         string; combining name of the level-1 folder and the condition_method
%       stats.channelpool:      string; summarizing the specified electrodes
%       stats.reduce_dims:      [] (only relevant for MVPA)
%       stats.cfg:              struct; the cfg of the input
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Common use cases:
%   (1) get ERP difference between conditions
%   (2) get condition-average ERP and test against baseline 
%   (3) get condition-specific ERPs and test each against baseline
%   (4) get N2pc/CDA electrode subtractions for two conditions, average them, and test against baseline
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
% cfg.electrode_method  = 'subtract'; --> specifies subtraction method, such that for condition 1 PO8 
%                                        is subtracted from P07, while for condition 2 P07 is
%                                        subtracted from PO8
%
% part of the ADAM toolbox, by J.J.Fahrenfort, VU, 2017/2018
% 
% See also ADAM_COMPUTE_GROUP_MVPA, ADAM_MVPA_FIRSTLEVEL, ADAM_PLOT_BDM_WEIGHTS

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

% Main routine, is a folder name specified? If not, pop up selection dialog
if isempty(folder_name)
    if ~isfield(cfg,'startdir')
        cfg.startdir = '';
        disp('NOTE: it is easier to select a directory when you indicate a starting directory using cfg.startdir, otherwise you have to start selection from root every time...');
    end
    folder_name = uigetdir(cfg.startdir,'select directory to plot');
    if ~ischar(folder_name)
        error('no folder was selected');
    end
    cfg.folder = folder_name;
    % where am I?
    ndirs = drill2data(folder_name);
    if isempty(plot_order)
        dirz = dir(folder_name);
        dirz = {dirz([dirz(:).isdir]).name};
        plot_order = dirz(cellfun(@isempty,strfind(dirz,'.')));
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
        dirz = dirz(cellfun(@isempty,strfind(dirz,'.')));
        if ~all(ismember(plot_order,dirz))
            error('One or more of the folders specified in cfg.plot_order cannot be found in this results directory. Change cfg.plot_order or select a different directory for plotting.');
        end
    end
    % loop through directories (results folders)
    stats = [];
    for cdirz = 1:numel(plot_order)
        stats = [stats subcompute_group_ERP(cfg,[folder_name filesep plot_order{cdirz}])];
    end
else
    if ~iscell(folder_name)
        folder_name = {folder_name};
    end
    stats = [];
    for cdirz=1:numel(folder_name)
        if ~exist(folder_name{cdirz},'dir')
            error([folder_name{cdirz} ' should refer to a full and existing folder path. Alternatively leave folder_name empty to pop up a selection dialog.']);
        end
        [stats, cfg] = [stats subcompute_group_ERP(cfg,folder_name{cdirz})];
    end
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
avg_conditions = false;
name = [];
reduce_dims = [];
mpcompcor_method = 'uncorrected';
timelim = [];
resample_eeg = 0;
electrode_def = [];
condition_def = [1,2]; % By default substracting cond1 - cond2
electrode_method = 'average';
condition_method = 'keep'; 
% unpack graphsettings
plottype = '2D';
channelpool = '';
v2struct(cfg);

% pack graphsettings with defaults
nameOfStruct2Update = 'cfg';
cfg = v2struct(tail,reduce_dims,indiv_pval,cluster_pval,name,plottype,mpcompcor_method,electrode_def,electrode_method,condition_def,condition_method,timelim,resample_eeg,nameOfStruct2Update);

% fill some empties
if isempty(electrode_def)
    error('no electrode_def was specified in cfg, set cfg.electrode_def to some electrode, e.g. ''Oz'', see help of this function for more info.');
end
pval(1) = indiv_pval;
pval(2) = cluster_pval;
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
cfg.plotFreq = plotFreq;
subjectfiles = dir([folder_name filesep channelpool plotFreq filesep '*.mat']);
[~, condname] = fileparts(folder_name);
subjectfiles = { subjectfiles(:).name };
nSubj = numel(subjectfiles);
if nSubj == 0
    error(['cannot find data in specified folder ' folder_name filesep channelpool plotFreq{:} ' maybe you should specify (a different) cfg.channelpool?']);
end

% prepare figure in case individual subjects are plotted
if plotsubjects;
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
    
    % compute electrodes and do electrode subtractions
    FT_ERP = restrict_FT_ERP(FT_ERP,cfg);
    settings.times = {FT_ERP.time, FT_ERP.time};
    
    % possible use cases:
    % 	(1) get ERP difference between conditions + stats against zero for one or more results folders
    % 	(2) get raw ERPs for one or more conditions from one or more results folders and average them
    %   (3) get raw ERPs for more conditions from a given results folder without averaging
    % 	(4) get N2pc/CDA electrode subtractions for two conditions, average
    %       them, and do stats against zero (for each results folder)
        
    % strmpi(electrode_method,'average')
    % strmpi(condition_method,'subtract')

    % now do condition subtraction
    clear trial;
    if strcmpi(condition_method,'subtract')
        if ~(size(condition_def,2)==2)
            error('Condition_def does not contain the correct number of conditions (2) to be able to subtract.');
        end
        trial = FT_ERP.trial(FT_ERP.trialinfo==condition_def(1),:,:) - FT_ERP.trial(FT_ERP.trialinfo==condition_def(2),:,:);
    elseif strcmpi(condition_method,'average')
        trial = mean(FT_ERP.trial(ismember(FT_ERP.trialinfo,condition_def),:,:),1);
    else  % not subtracting conditions, plain extraction
        trial = FT_ERP.trial(ismember(FT_ERP.trialinfo,condition_def),:,:);
    end
    for cCond=1:size(trial,1)
        if  all(all(all(trial<10^-3)))
            trial = trial*10^6; % little hack to change V to muV (if applicable)
        end
        ClassTotal{cCond}(cSubj,:) = squeeze(trial(cCond,:));
    end
    
    % plot individual subjects
    if plotsubjects
        subplot(numSubplots(nSubj,1),numSubplots(nSubj,2),cSubj);
        onestat.ClassOverTime = squeeze(trial);
        onestat.StdError = [];
        onestat.pVals = zeros(size(squeeze(trial)));
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
        tmpcfg.plot_model = plot_model;
        tmpcfg.plotsubjects = true;
        tmpcfg.plot_order = {condname};
        tmpcfg.acclim2D = [];
        tmpcfg.acclim3D = [];
        tmpcfg.acctick = [];
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
    if strcmpi(electrode_method,'subtract')
        condname =  [origcondname ' channel subtraction'];
    elseif strcmpi(condition_method,'subtract')
        condname =  [origcondname ' subtraction'];
    elseif strcmpi(condition_method,'average')
        condname = [origcondname ' average'];
    else
        condname = [origcondname ' erp' num2str(condition_def(cCond))];
    end
    
    % get some stats
    ClassAverage = mean(ClassTotal{cCond},1);
    ClassStdErr = std(ClassTotal{cCond},0,1)/sqrt(size(ClassTotal{cCond},1));

    if strcmp(mpcompcor_method,'fdr')
        % FDR CORRECTION
        [~,ClassPvals] = ttest(ClassTotal{cCond},chance,'tail',tail);
        thresh = fdr(squeeze(ClassPvals{cCond}),pval(2));
        ClassPvals(ClassPvals>thresh) = 1;
    elseif strcmp(mpcompcor_method,'cluster_based')
        % CLUSTER BASED CORRECTION
        [ClassPvals, pStruct] = cluster_based_permutation(ClassTotal{cCond},chance,cfg,settings);
        % compute Pstruct
    elseif strcmp(mpcompcor_method,'uncorrected')
        % NO MP CORRECTION
        [~,ClassPvals] = ttest(ClassTotal{cCond},chance,pval(1),tail);
    else
        % NO TESTING, PLOT ALL
        ClassPvals = zeros(1,size(ClassTotal{cCond},2));
    end
    
    % outputs: put it in a matrix for consistency in plot function
    settings.measuremethod = '\muV';
    stats(cCond).ClassOverTime = ClassAverage;
    stats(cCond).StdError = ClassStdErr;
    stats(cCond).pVals = ClassPvals;
    stats(cCond).mpcompcor_method = mpcompcor_method;
    stats(cCond).settings = settings;
    stats(cCond).condname = condname;
    stats(cCond).channelpool = FT_ERP.channelpool;
    if exist('pStruct','var')
        stats(cCond).pStruct = pStruct;
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
% subtracting electrode sets, subtracts electrode 2 from electrode 1 for each condition
if strcmpi(electrode_method,'subtract') %iscell(electrode_def{1}) 
    if size(electrode_def,2)~=2
        error('to subtract electrode sets, you should define two columns in cfg.electrode_def');
    end
    if size(electrode_def,1) ~= numel(condition_def)
        electrode_def = repmat(electrode_def,[numel(condition_def),1]);
    end
    for cCond = 1:size(electrode_def,1)
        for cDif=1:size(electrode_def,2)
            cfg = [];
            FT_TEMP(cDif) = select_channels_from_FT_EEG(FT_EEG,electrode_def{cCond,cDif});
        end
        trial(cCond,:,:) = FT_TEMP(1).trial(FT_TEMP(1).trialinfo==condition_def(cCond),:,:) - FT_TEMP(2).trial(FT_TEMP(2).trialinfo==condition_def(cCond),:,:);
        channelpool{cCond} = [FT_TEMP(1).label '-' FT_TEMP(2).label];
    end
    FT_EEG.trial = trial;
    FT_EEG.channelpool = channelpool;
elseif strcmpi(electrode_method,'average') % extract and average
    if ~iscell(electrode_def)
        electrode_def = {electrode_def};
    end
    FT_EEG = select_channels_from_FT_EEG(FT_EEG,electrode_def);
    FT_EEG.channelpool = FT_EEG.label;
else
    error('Specify cfg.electrode_method as ''average'' (default) or ''subtract''');
end
if any(size(FT_EEG.trial)==0)
    dims = regexp(FT_EEG.dimord,'_','split');
    error(['There are no dimensions left in these fields, change cfg selection parameters for ''' cellarray2csvstring(dims(logical(size(FT_EEG.trial)==0))) '''.']);
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

