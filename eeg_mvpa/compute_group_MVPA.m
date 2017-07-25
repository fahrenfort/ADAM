function [stats,cfg] = compute_group_MVPA(folder_name,cfg,mask)
% function [stats,cfg] = compute_group_MVPA(folder_name,cfg,mask)
%
% Extracts group classification data, classifier weights, forward model
% parameters etc. Also does basic stats on the extracted conditions.
% Use this as input for plot functions such as plot_MVPA,
% plot_MVPA_weights, plot_CTF and plot_FEM_weights.
% folder_name is the folder that contains the condition data.
% If folder_name is left empty, a selection dialog will pop up for you to
% indicate the root folder in which the folders that contain the condition
% folders are located (use this if you want to compute more than one
% condition at once). The conditions will be plotted in subplots by
% plot_MVPA and plot_CTF cfg is a struct that contains some input
% settings for compute_group_MVPA and ensueing plot functions:
% cfg.one_two_tailed = 'two' (can also be 'one')
% cfg.indiv_pval = .05;
% cfg.cluster_pval = .05;
% cfg.mpcompcor_method = 'uncorrected' (default, can also be 'cluster_based', 'fdr' or 'none')
% cfg.plot_dim = 'time_time' or 'freq_time' (default: 'time_time')
% cfg.reduce_dims = 'diag', 'avtrain', 'avtest' or 'avfreq' or []
% cfg.trainlim = is the time limits over which to constrain the
% training data, in ms.
% cfg.testlim = the time limits over which to constrain the testing
% data, in ms.
% cfg.timelim constrains trainlim and testlim at once (takes precedence
% over trainlim and testlim
% cfg.freqlim = is the frequency limits over which to constrain the
% frequency dimension.
% cfg.channelpool = 'ALL' (can also be e.g. 'OCCIP', see
% classify_RAW_eeglab_data.m and classify_TFR_from_eeglab_data.m for more 
% options.
% cfg.exclsubj = index numbers of subjects to skip (not including
% these subjects in the results structs). No subjects are excluded when
% left empty (default).
% Specify whether to extract the BDM or the FEM data using
% cfg.plot_model = 'BDM' (default) or 'FEM'
% You can also (optionally) specify which conditions you want to extract
% (and in which order) using cfg.plot_order, as a cell array.
%
% Example: 
% cfg.timelim = [0.2 1.2];
% cfg.mpcompcor_method = 'cluster_based';
% cfg.startdir = '/Volumes/backup/WM_debunk_EEG';
% cfg.plot_order = { 'LOCATION_TASK','SEARCH_TASK'};
% cfg.plot_model = 'FEM';
% [stats, cfg] = compute_group_MVPA('',cfg);
% 
% By J.J.Fahrenfort, VU 2015, 2016, 2017

% First get some settings
if nargin<3
    mask = [];
end

plot_order = {};

% backwards compatibility
plot_dim = 'time_time'; % default, 'time_time' or 'freq_time'
v2struct(cfg);
if exist('plotmodel','var')
    plot_model = plotmodel;
    cfg.plot_model = plot_model;
    cfg = rmfield(cfg,'plotmodel');
end
if exist('get_dim','var')
    plot_dim = get_dim;
    cfg.plot_dim = plot_dim;
    cfg = rmfield(cfg,'get_dim');
end
if strcmpi(plot_dim,'frequency_time') || strcmpi(plot_dim,'time_frequency') || strcmpi(plot_dim,'time_freq')
    plot_dim = 'freq_time';
    cfg.plot_dim = plot_dim;
end
if exist('plotorder','var')
    plot_order = plotorder;
    cfg.plot_order = plot_order;
    cfg = rmfield(cfg,'plotorder');
end

% check freqlimits
if (strcmpi(plot_dim,'freq_time') && ~strcmpi(reduce_dims,'avfreq')) && (numel(freqlim) == 1 || abs(diff(freqlim)) <= 2)
    wraptext('WARNING: your cfg.freqlim indicates a rather small range of frequencies given cfg.plot_dim ''freq_time'', use cfg.reduce_dims = ''avfreq'' if you intend to average. Now simply plotting all frequencies.');
    freqlim = [];
    cfg.freqlim = [];
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
%     if ~exist('plot_order','var') || isempty(plot_order)
%         dirz = dir(folder_name);
%         dirz = {dirz([dirz(:).isdir]).name};
%         dirz = dirz(cellfun(@isempty,strfind(dirz,'.'))); 
%     else
%         dirz = cfg.plot_order;
%     end
%     % loop through directories
%     for cdirz = 1:numel(dirz)
%         [stats(cdirz), cfg] = subcompute_group_MVPA([folder_name filesep dirz{cdirz}],cfg,mask);
%     end
    if isempty(plot_order)
        dirz = dir(folder_name);
        dirz = {dirz([dirz(:).isdir]).name};
        plot_order = dirz(cellfun(@isempty,strfind(dirz,'.')));
        % determine whether to drill down or not
        ndirs = drill2data(folder_name);
        if ndirs == 1
            [folder_name, plot_order] = fileparts(folder_name);
            plot_order = {plot_order};            
        elseif ndirs > 2
            error('You seem to be selecting a directory that is too high in the hiearchy, drill down a little more.');
        end
        cfg.plot_order = plot_order;
    end
    % loop through directories (results folders)
    for cdirz = 1:numel(plot_order)
        if numel(plot_order) == 1 % getting from single folder
            [stats, cfg] = subcompute_group_MVPA([folder_name filesep plot_order{cdirz}], cfg, mask);
        else % getting from multiple folders
            [stats(cdirz), cfg] = subcompute_group_MVPA([folder_name filesep plot_order{cdirz}], cfg, mask);
        end
    end
    cfg.folder = folder_name;
else
    if ~exist('folder_name','dir') && ~iscell(folder_name) 
        error([folder_name ' should refer to a full and existing folder path. Alternatively leave folder_name empty to pop up a selection dialog.']);
    end
    [stats, cfg] = subcompute_group_MVPA(folder_name,cfg,mask);
end

% subroutine for each condition
function [stats,cfg] = subcompute_group_MVPA(folder_name,cfg,mask)
% set defaults
one_two_tailed = 'two';
indiv_pval = .05;
cluster_pval = .05;
plotsubjects = false;
name = [];
channelpool = 'ALL'; % 'ALL', 'OCCIP', 'PARIET' etc, see select_channels.m function to make adjustments
mpcompcor_method = 'uncorrected';
plot_model = 'BDM'; % 'BDM' or 'FEM'
reduce_dims = []; % 'diag' 'avtrain' 'avtest' or 'avfreq'
timelim = [];
trainlim = [];
testlim = [];
freqlim = [];
exclsubj= [];
v2struct(cfg);
% general time limit
if ~isempty(timelim) % timelim takes precedence
    trainlim = timelim;
    testlim = timelim;
    cfg.trainlim = trainlim;
    cfg.testlim = testlim;
end

% unpack settings
v2struct(cfg);

% set defaults
pval(1) = indiv_pval;
pval(2) = cluster_pval;

% some logical checking: is this a frequency folder?
freqfolder_contains_time_time = ~isempty(dir([folder_name filesep channelpool filesep 'freq*'])); 
freqfolder_contains_freq_time = ~isempty(dir([folder_name filesep channelpool filesep 'allfreqs']));
freqfolder = any([freqfolder_contains_time_time freqfolder_contains_freq_time]);
if freqfolder
    if strcmpi(plot_dim,'freq_time') && freqfolder_contains_freq_time
        plotFreq = {[filesep 'allfreqs']};
    elseif strcmpi(plot_dim,'freq_time')
        disp('WARNING: freq_time is not availaible in this folder, defaulting to cfg.plot_dim = ''time_time''');
        plot_dim = 'time_time';
    end
    if strcmpi(plot_dim,'time_time') && freqfolder_contains_time_time
        if isempty(freqlim)
            freqlim = input('What frequency or frequency range should I extract (e.g. type [2 10] to average between 2 and 10 Hz)? ');
        end
        if numel(freqlim) > 1 % make a list of frequencies over which to average
            freqlist = dir([folder_name filesep channelpool filesep 'freq*']); freqlist = {freqlist(:).name};
            for c=1:numel(freqlist); freqsindir(c) = string2double(regexprep(freqlist{c},'freq','')); end
            freqsindir = sort(freqsindir); freqs2keep = find(freqsindir >= min(freqlim) & freqsindir <= max(freqlim));
            for c=1:numel(freqs2keep); plotFreq{c} = [filesep 'freq' num2str(freqsindir(freqs2keep(c)))]; end
        else % or just a single frequency
            plotFreq{1} = [filesep 'freq' num2str(freqlim)];
        end
    elseif strcmpi(plot_dim,'time_time') && ~freqfolder_contains_time_time
        disp('WARNING: time_time is not available in this folder, defaulting to cfg.plot_dim = ''freq_time''');
        plot_dim = 'freq_time';
        plotFreq{1} = [filesep 'allfreqs'];
    end
else
    plotFreq = {''};
    freqlim = [];
end
        
% pack graphsettings with defaults
nameOfStruct2Update = 'cfg';
cfg = v2struct(freqlim,plotFreq,trainlim,testlim,one_two_tailed,indiv_pval,cluster_pval,plot_model,mpcompcor_method,reduce_dims,freqlim,nameOfStruct2Update);

% get filenames
subjectfiles = dir([folder_name filesep channelpool plotFreq{1} filesep '*.mat']);
[~, condname] = fileparts(folder_name);
subjectfiles = { subjectfiles(:).name };

% limiting subjects
if ~isempty(exclsubj)
   subjectfiles = select_subjects(subjectfiles,exclsubj,true);
end

% see if data exists
nSubj = numel(subjectfiles);
if nSubj == 0
    error(['cannot find data in specified folder ' folder_name filesep channelpool plotFreq]);
end

% prepare figure in case individual subjects are plotted
if plotsubjects;
    fh = figure('name',['individual subjects, condition: ' condname]);
    set(fh, 'Position', get(0,'Screensize'));
    set(fh,'color','w');
end

% do the loop, restrict time and frequency if applicable
firstchanlocs = [];
for cSubj = 1:nSubj
    fprintf(1,'loading subject %d of %d\n', cSubj, nSubj);
    
    % initialize subject
    clear ClassOverTimeAv WeightsOverTimeAv covPatternsOverTimeAv corPatternsOverTimeAv C2_averageAv C2_perconditionAv;

    % loop over frequencies (if no frequencies exist, it simply loads raw)
    for cFreq = 1:numel(plotFreq)
        
        % locate data
        matObj = matfile([folder_name filesep channelpool plotFreq{cFreq} filesep subjectfiles{cSubj}]);
        settings = matObj.settings;
        
        % for backward compatibility
        if strcmpi(settings.dimord,'frequency_time')
            settings.dimord = 'freq_time';
        end
        v2struct(settings);
        
        % get data
        if ~isempty(whos(matObj,'BDM')) && strcmpi(plot_model,'BDM')
            v2struct(matObj.BDM); % unpack 
        elseif ~isempty(whos(matObj,'FEM')) && strcmpi(plot_model,'FEM')
            v2struct(matObj.FEM); % unpack 
        else
            error('cannot find data');
        end
        
        % find limits
        if ~exist('firstchanlocs','var')
            firstchanlocs = [];
        end
        [settings, cfg, lim1, lim2, dataindex, firstchanlocs] = find_limits(settings, cfg, firstchanlocs);
        v2struct(cfg);
        
        % limit ClassOverTime
        ClassOverTime = ClassOverTime(lim1,lim2);
        
        % limit weights too
        if strcmpi(dimord,'freq_time')
            WeightsOverTime = WeightsOverTime(lim1,lim2,dataindex,:);
        else
            WeightsOverTime = WeightsOverTime(lim2,dataindex,:);
        end
        if strcmpi(plot_model,'BDM')
            if strcmpi(dimord,'freq_time')
                covPatternsOverTime = covPatternsOverTime(lim1,lim2,dataindex);
                corPatternsOverTime = corPatternsOverTime(lim1,lim2,dataindex);
            else
                covPatternsOverTime = covPatternsOverTime(lim2,dataindex);
                corPatternsOverTime = corPatternsOverTime(lim2,dataindex);
            end
        else
            C2_average = C2_average(lim1,lim2,:);
            C2_percondition = C2_percondition(lim1,lim2,:,:);
        end
        
        % if applicable, reduce dimensionality (creates 2D plot)
        if strcmpi(reduce_dims,'avfreq') && strcmpi(dimord,'freq_time')
            if isempty(freqlim)
                disp('WARNING: you are averaging across ALL frequencies, are you sure that is what you want?');
            end
            ClassOverTime = mean(ClassOverTime,1);
            WeightsOverTime = mean(WeightsOverTime,1);
            if strcmpi(plot_model,'BDM')
                covPatternsOverTime = mean(covPatternsOverTime,1);
                corPatternsOverTime = mean(corPatternsOverTime,1);
%           NOTE: the plot_CTF function still assumes the full matrix, 
%           needs to be updated. For now just pass the full matrix. When
%           this is fixed, could also plot CTF across al testing points
%           when training on one specific timepoint or vice versa
%             else
%                 C2_average = mean(C2_average,1);
%                 C2_percondition = mean(C2_percondition,1);
            end
            mask = sum(mask,1);
        elseif strcmpi(reduce_dims,'avtrain') && strcmpi(dimord,'time_time')
            if isempty(trainlim)
                disp('WARNING: you are averaging across ALL training time points, are you sure that is what you want?');
            end
            ClassOverTime = mean(ClassOverTime,2); % IMPORTANT, TRAIN IS ON SECOND DIMENSION
%             if strcmpi(plot_model,'FEM')
%                 C2_average = mean(C2_average,2);
%                 C2_percondition = mean(C2_percondition,2);
%             end
            mask = sum(mask,2);
        elseif strcmpi(reduce_dims,'avtest') && strcmpi(dimord,'time_time')
            if isempty(trainlim)
                disp('WARNING: you are averaging across ALL testing time points, are you sure that is what you want?');
            end
            ClassOverTime = mean(ClassOverTime,1); % IMPORTANT, TEST IS ON FIRST DIMENSION
%             if strcmpi(plot_model,'FEM')
%                 C2_average = mean(C2_average,1);
%                 C2_percondition = mean(C2_percondition,1);
%             end
            mask = sum(mask,1);
        elseif strcmpi(reduce_dims,'diag') && strcmpi(dimord,'time_time')
            ClassOverTime = diag(ClassOverTime);
%             if strcmpi(plot_model,'FEM')
%                 for c1 =1:size(C2_percondition,3)
%                     diagC2_average(:,c1) = diag(C2_average(:,:,c1));
%                     for c2 = 1:size(C2_percondition,4)
%                         diagC2_percondition(:,c1,c2) = diag(C2_percondition(:,:,c1,c2));
%                     end
%                 end
%                 C2_average = diagC2_average;
%                 C2_percondition = diagC2_percondition;
%             end
            mask = diag(mask);
        elseif strcmpi(reduce_dims,'diag') && strcmpi(dimord,'freq_time')
            disp('WARNING: cannot reduce dimensionality along diagonal when dimord is freq_time');
        end
        
        % sum up to compute average over frequencies (avfreq)
        if ~exist('ClassOverTimeAv','var'); ClassOverTimeAv = zeros(size(ClassOverTime)); end
        if ~exist('WeightsOverTimeAv','var'); WeightsOverTimeAv = zeros(size(WeightsOverTime)); end
        ClassOverTimeAv = ClassOverTimeAv + ClassOverTime;
        WeightsOverTimeAv = WeightsOverTimeAv + WeightsOverTime;
        if strcmpi(plot_model,'BDM')
            if ~exist('covPatternsOverTimeAv','var'); covPatternsOverTimeAv = zeros(size(covPatternsOverTime)); end
            if ~exist('corPatternsOverTimeAv','var'); corPatternsOverTimeAv = zeros(size(corPatternsOverTime)); end
            covPatternsOverTimeAv = covPatternsOverTimeAv + covPatternsOverTime;
            corPatternsOverTimeAv = corPatternsOverTimeAv + corPatternsOverTime;
        else
            if ~exist('C2_averageAv','var'); C2_averageAv = zeros(size(C2_average)); end
            if ~exist('C2_perconditionAv','var'); C2_perconditionAv = zeros(size(C2_percondition)); end
            C2_averageAv = C2_averageAv + C2_average;
            C2_perconditionAv = C2_perconditionAv + C2_percondition;
        end
        
    end
    
    % by default it computes the average over frequencies when specifying
    % time_time (cfg.reduce_dims = 'avfreq' is actually superfluous in this case)
    ClassOverTime = ClassOverTimeAv / numel(plotFreq);
    WeightsOverTime = WeightsOverTimeAv / numel(plotFreq);
    if strcmpi(plot_model,'BDM')
        covPatternsOverTime = covPatternsOverTimeAv / numel(plotFreq);
        corPatternsOverTime = corPatternsOverTimeAv / numel(plotFreq);
    else
        C2_average = C2_averageAv / numel(plotFreq);
        C2_percondition = C2_perconditionAv / numel(plotFreq);
    end
    
    % make big matrix of of all subjects
    ClassOverTimeAll{1}(cSubj,:,:) = ClassOverTime;
    indx = [{cSubj} repmat({':'}, 1, ndims(WeightsOverTime))];
    WeightsOverTimeAll(indx{:}) = WeightsOverTime;
    if strcmpi(plot_model,'BDM')
        indx = [{cSubj} repmat({':'}, 1, ndims(covPatternsOverTime))];
        covPatternsOverTimeAll(indx{:}) = covPatternsOverTime;
        indx = [{cSubj} repmat({':'}, 1, ndims(corPatternsOverTime))];
        corPatternsOverTimeAll(indx{:}) = corPatternsOverTime;
    else
        indx = [{cSubj} repmat({':'}, 1, ndims(C2_average))];
        C2_averageAll(indx{:}) = C2_average;
        indx = [{cSubj} repmat({':'}, 1, ndims(C2_percondition))];
        C2_perconditionAll(indx{:}) = C2_percondition;
    end
    
    % plot individual subjects
    if plotsubjects
        subplot(numSubplots(nSubj,1),numSubplots(nSubj,2),cSubj);
        onestat.ClassOverTime = ClassOverTime;
        onestat.StdError = [];
        onestat.pVals = zeros(size(ClassOverTime));
        onestat.indivClassOverTime = [];
        onestat.settings = settings;
        onestat.condname = condname;
        onestat.channelpool = channelpool;
        tmpcfg = cfg;
        tmpcfg.plotsubject = true;
        tmpcfg.plot_order = {condname};
        plot_MVPA(onestat,tmpcfg);
        subjname = subjectfiles{cSubj};
        underscores = strfind(subjname,'_');
        subjname = regexprep(subjname(underscores(2)+1:underscores(end)-1),'_',' ');
        ntitle(subjname,'fontsize',10,'fontweight','bold');
    end
end

% determine chance level
if strcmpi(settings.measuremethod,'hr-far') || strcmpi(plot_model,'FEM')
    chance = 0;
else
    chance = 1/settings.nconds;
end

% compute standard errors and averages
ClassStdErr(1:size(ClassOverTimeAll{1},2),1:size(ClassOverTimeAll{1},3)) = std(ClassOverTimeAll{1},0,1)/sqrt(size(ClassOverTimeAll{1},1));
if sum(sum(ClassStdErr)) == 0 ClassStdErr = []; end % don't plot stderror when there is none
ClassAverage(1:size(ClassOverTimeAll{1},2),1:size(ClassOverTimeAll{1},3)) = mean(ClassOverTimeAll{1},1);
ClassOverTimeAll{2} = repmat(chance,size(ClassOverTimeAll{1}));

% statistical testing
if nSubj > 1
    if strcmpi(mpcompcor_method,'fdr')
        % FDR CORRECTION
        if strcmpi(one_two_tailed,'two')
            [~,ClassPvals(1:size(ClassOverTimeAll{1},2),1:size(ClassOverTimeAll{1},3))] = ttest(ClassOverTimeAll{1},ClassOverTimeAll{2},pval(1),'both');
        else
            [~,ClassPvals(1:size(ClassOverTimeAll{1},2),1:size(ClassOverTimeAll{1},3))] = ttest(ClassOverTimeAll{1},ClassOverTimeAll{2},pval(1),'right');
        end
        thresh = fdr(squeeze(ClassPvals),pval(2));
        ClassPvals(ClassPvals>thresh) = 1;
    elseif strcmpi(mpcompcor_method,'cluster_based')
        % CLUSTER BASED CORRECTION
        [ClassPvals, pStruct] = cluster_based_permutation(ClassOverTimeAll{1},ClassOverTimeAll{2},cfg,settings,mask);
    elseif strcmpi(mpcompcor_method,'uncorrected')
        % NO MP CORRECTION
        if strcmpi(one_two_tailed,'two')
            [~,ClassPvals(1:size(ClassOverTimeAll{1},2),1:size(ClassOverTimeAll{1},3))] = ttest(ClassOverTimeAll{1},ClassOverTimeAll{2},pval(1),'both');
        else
            [~,ClassPvals(1:size(ClassOverTimeAll{1},2),1:size(ClassOverTimeAll{1},3))] = ttest(ClassOverTimeAll{1},ClassOverTimeAll{2},pval(1),'right');
        end
        ClassPvals(~mask) = 1;
    else
        % NO TESTING, PLOT ALL
        ClassPvals = zeros([size(ClassOverTimeAll{1},2) size(ClassOverTimeAll{1},3)]);
    end
else
    ClassPvals = zeros([size(ClassOverTimeAll{1},2) size(ClassOverTimeAll{1},3)]);
end

% outputs
stats.ClassOverTime = ClassAverage;
stats.StdError = ClassStdErr;
stats.pVals = ClassPvals;
stats.mpcompcor_method = mpcompcor_method;
stats.indivClassOverTime = ClassOverTimeAll{1};
stats.settings = settings;
stats.condname = condname;
stats.filenames = subjectfiles;
stats.channelpool = channelpool;
if exist('pStruct','var')
    stats.pStruct = pStruct;
end
%cfg = v2struct(name,nameOfStruct2Update);

% compute weights stuff
if exist('WeightsOverTimeAll','var')
    weights.avWeights = squeeze(mean(WeightsOverTimeAll,1));
    weights.indivWeights = squeeze(WeightsOverTimeAll);
end
if strcmpi(plot_model,'BDM')
    weights.avCovPatterns = squeeze(mean(covPatternsOverTimeAll,1));   
    weights.indivCovPatterns = squeeze(covPatternsOverTimeAll);
    weights.avCorPatterns = squeeze(mean(corPatternsOverTimeAll,1));
    weights.indivCorPatterns = squeeze(corPatternsOverTimeAll);
else
    weights.CTF = squeeze(mean(C2_averageAll,1));
    weights.semCTF = squeeze(std(C2_averageAll,0,1)/sqrt(size(C2_averageAll,1)));
    weights.indivCTF = squeeze(C2_averageAll);
    CTFpercond = squeeze(mean(C2_perconditionAll,1));
    semCTFpercond = squeeze(std(C2_perconditionAll,0,1)/sqrt(size(C2_perconditionAll,1)));
    indivCTFpercond = squeeze(C2_perconditionAll);
    % make a nicer list of CTFs per condition, better for plotting later on
    for cCond = 1:size(CTFpercond,3)
        weights.CTFpercond{cCond} = squeeze(CTFpercond(:,:,cCond,:));
        weights.semCTFpercond{cCond} = squeeze(semCTFpercond(:,:,cCond,:));
        weights.indivCTFpercond{cCond} = squeeze(indivCTFpercond(:,:,:,cCond,:));
    end
end
stats.weights = weights;
stats.cfg = cfg;
disp('done!');

function [settings, cfg, lim1, lim2, dataindex, firstchanlocs] = find_limits(settings, cfg, firstchanlocs) 
% find limits within which to constrain ClassOverTime
v2struct(cfg); % unpack cfg
v2struct(settings); % unpack settings
if strcmpi(dimord,'freq_time') 
    freqs = settings.freqs; % to fix that freqs is also a function and v2struct can't deal with that
end

if numel(settings.times) == 1 && strcmpi(settings.dimord,'time_time')
    times{2} = times{1};
end

% get the relevant electrodes and obtain the correct order for weights
if ~isfield(settings,'chanlocs') % if no chanlocdata exist in settings
    if ~exist('chanlocdata','var')
        chanlocdata = readlocs('plotting_1005.sfp','importmode','native');
    end
    [~, chanindex, dataindex] = intersect({chanlocdata(:).labels},settings.channels,'stable');
    chanlocs = chanlocdata(chanindex); % put all in the same order as imported locations
else % otherwise just extract from settings
    chanlocs = settings.chanlocs;
    if iscell(chanlocs)
        chanlocs = chanlocs{1};
    end
    if isempty(firstchanlocs)
        firstchanlocs = chanlocs;
    end
    [~, ~, dataindex] = intersect({firstchanlocs(:).labels},{chanlocs(:).labels},'stable');
    chanlocs = firstchanlocs;
end
if numel(chanlocs) < numel(settings.channels)
    error('could not find location info for all channels');
end
settings.chanlocs = chanlocs;

% continue limit operation
% NOTE: ClassOverTime has dimensions: test_time * train_time OR freq * time
% In settings, times{1} is always train and times{2} is always test, but 
% ClassOverTime(1,:) is the first element of test_time (1st dimension) and
% ClassOverTime(:,1) is the first element of train_time (2nd dimension)
if strcmpi(dimord,'freq_time') && numel(freqlim)>1
    lim1 = nearest(freqs,min(freqlim)):nearest(freqs,max(freqlim));
    freqs = freqs(lim1);
elseif strcmpi(dimord,'freq_time') && numel(freqlim) == 1
    lim1 = nearest(freqs,freqlim);
    freqs = freqs(lim1);
elseif strcmpi(dimord,'freq_time') && isempty(freqlim)
    lim1 = true(size(freqs));
elseif strcmpi(dimord,'time_time') && ~isempty(testlim)
    lim1 = nearest(times{2}*1000,testlim(1)):nearest(times{2}*1000,testlim(2));
    times{2} = times{2}(lim1); % that is why times{2}(lim1)!
else
    lim1 = true(size(times{2})); % that is why lim1 = true(size(times{2}))!
end

% the time dimension (always present)
if ~isempty(trainlim)
    lim2 = nearest(times{1}*1000,trainlim(1)):nearest(times{1}*1000,trainlim(2));
    times{1} = times{1}(lim2); % that is why times{1}(lim2)!
else
    lim2 = true(size(times{1}));
end

% if the diagonal is plotted in 2D, restriction should be matched
if strcmpi(reduce_dims,'diag') && strcmpi(dimord,'time_time')
    lim1 = lim2;
end

% consolidate
settings.times = times;
if strcmpi(dimord,'freq_time')
    settings.freqs = freqs;
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