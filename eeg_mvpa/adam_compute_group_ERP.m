function [stats,cfg] = adam_compute_group_ERP(folder_name,cfg)
% function [stats,cfg] = adam_compute_group_ERP(folder_name,cfg)
% Computes group ERPs and extracts ERP averages
%
% Common use cases:
% 	(1) get ERP difference between conditions + stats against zero for one or more results folders
% 	(2) get raw ERPs for one or more conditions from one or more results folders and average them
%   (3) get raw ERPs for more conditions from a given results folder without averaging
% 	(4) get N2pc/CDA electrode subtractions for two conditions, average
%       them, and do stats against zero (for each results folder) 
% 
% Use case (1):
% specify electrodes to be extracted and averaged:
% cfg.electrode_def = { 'Oz', 'Iz', 'POz' };
%   - a cell array with the electrodes which will be averaged 
% cfg.elecrode_method = 'average' (default)
%   - averages signal over electrodes
% cfg.condition_def = [1,2];
% cfg.condition_methods = 'subtract'
%   - condition_def = [1,2] will subtract 2 from 1 and test them against
%     each other -> can do this for multiple results folders at once
%
% Use case (2):
% cfg.electrode_def = { 'Oz' };
%   - a cell array with the electrodes which will be averaged (can also be
%     a single electrode)
% cfg.condition_def = [1,2,3,4];
% cfg.condition_methods = 'average'
%   - will output the average of all conditions and test against 0 
%   -> can do this for multiple results folders at once
%
% Use case (3):
% cfg.electrode_def = { 'Oz' };
%   - a cell array with the electrodes which will be averaged (can also be
%     a single electrode)
% cfg.condition_def = [1,2,3,4];
% cfg.condition_methods = 'keep' (default)
%   - will output those conditions
%   -> can only do this for a single results folders
%
% Use case (4):
% cfg.electrode_def = {{'PO7'},{'PO8'};{'PO8'},{'PO7'}};
%   - a cell array with cell arrays of electrodes to be subtracted
% cfg.electrode_method = 'subtract'
%   - specifies subtraction method, such that for condition 1 PO8 is
%   subtracted from P07, while for condition 2 P07 is subtracted from PO8
% cfg.condition_def = [1,2];
% cfg.condition_methods = 'average'
%   - will output a stats structure for the average tested against 0
%   -> can do this for multiple results folders at once
%
% cfg can also specify the time interval, and the correction method and threshold for statistics
%
% Use stats output from this function as input for plot function adam_plot_MVPA
%
% By J.J.Fahrenfort, VU, 2016, 2017

plot_order = {};

% backwards compatibility
v2struct(cfg);
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
    for cdirz = 1:numel(plot_order)
        if numel(plot_order) == 1 % getting from single folder
            [stats, cfg] = subcompute_group_ERP([folder_name filesep plot_order{cdirz}],cfg);
        else % getting from multiple folders
            [stats(cdirz), cfg] = subcompute_group_ERP([folder_name filesep plot_order{cdirz}],cfg);
        end
    end
else
    if ~exist('folder_name','dir') && ~iscell(folder_name) 
        error([folder_name ' should refer to a full and existing folder path. Alternatively leave folder_name empty to pop up a selection dialog.']);
    end
    [stats, cfg] = subcompute_group_ERP(folder_name,cfg);
end
% fill cfg.plot_order in case not given by user and from single folder
if numel(cfg.plot_order) == 1
    cfg.plot_order = {stats(:).condname};
end

function [stats,cfg] = subcompute_group_ERP(folder_name,cfg)
% set defaults
one_two_tailed = 'two';
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
condition_method = '';
% unpack graphsettings
plottype = '2D';
channelpool = '';
v2struct(cfg);

% pack graphsettings with defaults
nameOfStruct2Update = 'cfg';
cfg = v2struct(one_two_tailed,reduce_dims,indiv_pval,cluster_pval,name,plottype,mpcompcor_method,electrode_def,electrode_method,condition_def,condition_method,timelim,resample_eeg,nameOfStruct2Update);

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
        adam_plot_MVPA(onestat,tmpcfg);
        subjname = subjectfiles{cSubj};
        underscores = strfind(subjname,'_');
        subjname = regexprep(subjname(underscores(2)+1:underscores(end)-1),'_',' ');
        ntitle(subjname,'fontsize',10,'fontweight','bold');
        drawnow;
    end
end

% next, compute stats
chance = 0;
if strcmp(one_two_tailed,'two')
    tail = 'both';
else
    tail = 'right';
end

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
        condname = ['condition ' num2str(condition_def(cCond))];
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
    error('Specify cfg.electrode_method as ''average''  or ''subtract''');
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

