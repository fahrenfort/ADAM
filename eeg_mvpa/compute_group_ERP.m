function [stats,cfg] = compute_group_ERP(folder_name,cfg)
% function [stats,weights,cfg] = compute_group_ERP(folder_name,cfg)
% electrode_sets = {{'PO7'},{'PO8'}}; 
%   - a cell array with the electrodes to be extracted and averaged { 'Oz', 'Iz', 'POz' };
%   - a cell array of cell arrays for which the electrodes will be extracted, averaged and substracted
%     e.g. in definition {{'P8'},{'P7'}} substracted will be {{'P7'} - {'P8'}}
%     in definition {{'PO7'},{'PO8'};{'PO8'},{'PO7'}}; will subtract PO8
%     from PO7 in the 1st condition and PO7 from PO8 in the 2nd condition
% condition_def = [1,2; 3,4 ];
%   If only one column will extract those conditions, if
%   cfg.avg_conditions = true it will also average those conditions
%   (default: false)
%   If two columns, conditions will be substracted using column1 - column2
%   If condition_def is not specified it will substract condition 2 from
%   condition 1 (default).
%   e.g. condition_def = [1;2;3;4;5]; will extract condition 1:5
%   if cfg.avg_conditions = true it will average them to one time
%   series
%   if condition_def = [1,2;3,4] it will subtract 2 from 1 and 4 from 3
%   to create only two time series
% cfg specifies the time interval, also does statistics against 0
% Computes group ERPs and extracts ERP averages
% Use this as input for plot function plot_MVPA
%
% By J.J.Fahrenfort, VU, 2016

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
        return
    end
    if ~exist('plot_order','var') || isempty(plot_order)
        dirz = dir(folder_name);
        dirz = {dirz([dirz(:).isdir]).name};
        dirz = dirz(cellfun(@isempty,strfind(dirz,'.'))); 
    else
        dirz = cfg.plot_order;
    end
    % loop through directories
    for cdirz = 1:numel(dirz)
        if numel(cfg.plot_order) == 1 % getting from single folder
            [stats, cfg] = subcompute_group_ERP([folder_name filesep dirz{cdirz}],cfg);
        else % getting from multiple folders
            [stats(cdirz), cfg] = subcompute_group_ERP([folder_name filesep dirz{cdirz}],cfg);
        end
    end
    cfg.folder = folder_name;
else
    if ~exist('folder_name','dir') && ~iscell(folder_name) 
        error([folder_name ' should refer to a full and existing folder path. Alternatively leave folder_name empty to pop up a selection dialog.']);
    end
    [stats, cfg] = subcompute_group_ERP(folder_name,cfg);
end

% % Main routine, is a folder name specified? If not, pop up selection dialog
% if isempty(folder_name)
%     if ~isfield(cfg,'startdir')
%         cfg.startdir = '';
%         disp('NOTE: it is easier to run this function if you indicate a starting directory in cfg.startdir');
%     end
%     folder_name = uigetdir(cfg.startdir,'select directory to plot');
%     if ~ischar(folder_name)
%         stats = [];
%         weights = [];
%         return
%     end
%     if ~isfield(cfg,'plotorder') || isempty(cfg.plotorder)
%         dirz = dir(folder_name);
%         dirz = {dirz([dirz(:).isdir]).name};
%         dirz = dirz(cellfun(@isempty,strfind(dirz,'.')));
%         cfg.plotorder = dirz;
%     else
%         dirz = cfg.plotorder;
%     end
%     for cdirz = 1:numel(dirz)
%         if numel(cfg.plotorder) == 1 % getting from single folder
%             [stats, cfg] = subcompute_group_ERP([folder_name filesep dirz{cdirz}],cfg);
%         else % getting from multiple folders
%             [stats(cdirz), cfg] = subcompute_group_ERP([folder_name filesep dirz{cdirz}],cfg);
%         end
%     end
% else
%     if ~exist('folder_name','dir')
%         error([folder_name ' should refer to a full and existing folder path']);
%     end
%     [stats, cfg] = subcompute_group_ERP(folder_name,cfg);
% end
% 
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
electrode_sets = [];
condition_def = [1,2]; % By default substracting cond1 - cond2
% unpack graphsettings
v2struct(cfg);
plottype = '2D';
channelpool = 'ALL';
if size(condition_def,2) == 2 && size(electrode_sets,2) == 2
    average = true; % average conditions when computing average N2pcs (so comparing left and right, make sure you don't do contra - ipsi, but use same elec1 - elec2 subtraction for both hemispheres
else
    average = false;
end

% pack graphsettings with defaults
nameOfStruct2Update = 'cfg';
cfg = v2struct(one_two_tailed,reduce_dims,indiv_pval,cluster_pval,name,plottype,mpcompcor_method,timelim,electrode_sets,resample_eeg,condition_def,nameOfStruct2Update);

% fill some empties
if isempty(electrode_sets)
    error('no electrode_sets was specified in settings, set cfg.electrode_sets to some electrode, e.g. ''Oz''.');
end
pval(1) = indiv_pval;
pval(2) = cluster_pval;

% get filenames
plotFreq = ''; % this is empty for now, but might be used to look at the ERPs obtained from a TF analysis
cfg.plotFreq = plotFreq;
subjectfiles = dir([folder_name filesep channelpool plotFreq filesep '*.mat']);
[~, condname] = fileparts(folder_name);
name = strsplit(folder_name,filesep);
name = name(1:end-1);
subjectfiles = { subjectfiles(:).name };
nSubj = numel(subjectfiles);
if nSubj == 0
    error(['cannot find data in specified folder ' folder_name filesep channelpool plotFreq filesep]);
end

% do the loop, restrict time and frequency if applicable
for cSubj = 1:nSubj
    fprintf(1,'loading subject %d of %d\n', cSubj, nSubj);
    matObj = matfile([folder_name filesep channelpool plotFreq filesep subjectfiles{cSubj}]);
    if ~isempty(whos(matObj,'settings'))
        settings = matObj.settings;
    else % BW compatibility, will become obsolete over time
        error('Cannot find settings struct which the contains the ERP info');
    end
    % OLD for backwards compatible
    if isfield(settings,'ft_erpstruct')
        FT_ERP = settings.ft_erpstruct;
    elseif isfield(settings,'FT_ERP')
        FT_ERP = settings.FT_ERP;
    else
        error('no ERP was computed for these data');
    end
    
    if iscell(FT_ERP)
        FT_ERP = FT_ERP{2};
        disp('Extracting ERPs from testing data');
    end
    
    % compute electrodes and do electrode substractions
    FT_ERP = restrict_FT_ERP(FT_ERP,cfg);
    settings.times = {FT_ERP.time, FT_ERP.time};
    
    % now do condition subtraction
    if size(condition_def,2) == 2 % subtracting
        for cCond=1:size(condition_def,1)
            if average == true
                trial(cCond,:,:) = squeeze(FT_ERP.trial(condition_def(cCond,1),:,:) - FT_ERP.trial(condition_def(cCond,2),:,:))/2;
            else
                trial(cCond,:,:) = squeeze(FT_ERP.trial(condition_def(cCond,1),:,:) - FT_ERP.trial(condition_def(cCond,2),:,:));
            end
        end
    else % not subtracting conditions, plain extraction
        if average % average over all conditions
            trial = mean(FT_ERP.trial(condition_def,:,:),1);
        else % per condition
            trial = squeeze(FT_ERP.trial(condition_def,:,:));
        end
    end
    for cCond=1:size(trial,1)
        ClassTotal{cCond}(cSubj,:) = squeeze(trial(cCond,:));
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
for cCond = 1:size(condition_def,1) % loop over conditions
    
    % determine cond_name when getting single folder (condition_def is given by user)
    if size(condition_def,2) == 2 && numel(plotorder) == 1 % subtracting
        condname = [ num2str(condition_def(cCond,1)) '-' num2str(condition_def(cCond,2))];
        if average
            channeldef = ['((' FT_ERP.channelpool{condition_def(cCond,1)} ') - (' FT_ERP.channelpool{condition_def(cCond,2)} ')) / 2'];
        else
            channeldef = ['(' FT_ERP.channelpool{condition_def(cCond,1)} ') - (' FT_ERP.channelpool{condition_def(cCond,2)} ')'];
        end
    %elseif size(condition_def,2) == 1 && numel(plotorder) == 1 % not subtracting
    %    condname = num2str(condition_def(cCond));
    %    channeldef = FT_ERP.channelpool{condition_def(cCond)};
    else
        channeldef = FT_ERP.channelpool; %FT_ERP.channelpool{condition_def(cCond)};
        condname = ['condition ' num2str(cCond)];
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
        [~,ClassPvals] = ttest(ClassTotal{cCond},chance,'tail',tail);
    else
        % NO TESTING, PLOT ALL
        ClassPvals = zeros(1,size(ClassTotal{cCond},2));
    end
    
    % outputs: put it in a matrix for consistency in plot function
    settings.measuremethod = '\muV';
    stats(cCond).ClassOverTime = ClassAverage;
    stats(cCond).StdError = ClassStdErr;
    stats(cCond).pVals = ClassPvals;
    stats(cCond).settings = settings;
    stats(cCond).condname = condname;
    stats(cCond).channelpool = channeldef;
    if exist('pStruct','var')
        stats(cCond).pStruct = pStruct;
    end
end
cfg = v2struct(name,nameOfStruct2Update);

function [FT_EEG] = restrict_FT_ERP(FT_EEG,cfg)
% resample / restrict the ERP
v2struct(cfg);
% resample?
if resample_eeg
    cfg = [];
    cfg.resamplefs = resample_eeg;
    FT_EEG = ft_resampledata(cfg,FT_EEG);
end
% limit time?
if ~isempty(timelim)
    cfg = [];
    cfg.latency = timelim/1000; % should be in seconds
    FT_EEG = ft_selectdata(cfg,FT_EEG);
end
clear channelpool;
% subtracting electrode sets, subtracts electrode 2 from electrode 1 for each condition
if iscell(electrode_sets{1}) 
    if size(electrode_sets,2)~=2
        error('to subtract electrode sets, you should define two columns in cfg.electrode_sets');
    end
    conds = unique(cfg.condition_def);
    if size(electrode_sets,1) ~= numel(conds)
        electrode_sets = repmat(electrode_sets,[numel(conds),1]);
    end
    for cCond = 1:size(electrode_sets,1)
        for cDif=1:size(electrode_sets,2)
            cfg = [];
            cfg.channel = electrode_sets{cCond,cDif};
            cfg.avgoverchan = 'yes';
            FT_TEMP{cDif} = ft_selectdata(cfg,FT_EEG);
        end
        trial(cCond,:,:) = FT_TEMP{1}.trial(cCond,:,:) - FT_TEMP{2}.trial(cCond,:,:);
        channelpool{cCond} = [electrode_sets{cCond,1}{:} '-' electrode_sets{cCond,2}{:}];
    end
    FT_EEG.trial = trial;
    FT_EEG.channelpool = channelpool;
else % plain extraction
    cfg.channel = electrode_sets;
    FT_EEG = ft_selectdata(cfg,FT_EEG);
    FT_EEG.channelpool = electrode_sets;
end

    
