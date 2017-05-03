function [stats,weights,gsettings] = compute_group_MVPA(folder_name,gsettings,mask)
% function [stats,weights,gsettings] = compute_group_MVPA(folder_name,gsettings,mask)
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
% plot_MVPA and plot_CTF gsettings is a struct that contains some input
% settings for compute_group_MVPA and ensueing plot functions:
% gsettings.one_two_tailed = 'two' (can also be 'one')
% gsettings.indiv_pval = .05;
% gsettings.cluster_pval = .05;
% gsettings.mpcompcor_method = 'uncorrected' (default, can also be 'cluster_based', 'fdr' or 'none')
% gsettings.plottype = '2D' (default, can also be '3D'), determines whether
% only to extract and test the diagonal, or whether to etract time_time or time_frequency
% gsettings.trainlim = is the time limits over which to constrain the
% training data, in ms.
% gsettings.testlim = the time limits over which to constrain the testing
% data, in ms.
% gsettings.timelim constrains trainlim and testlim at once (only works if
% trainlim and testlim are empty)
% gsettings.freqlim = is the frequency limits over which to constrain the
% frequency dimension.
% gsettings.channelpool = 'ALL' (can also be e.g. 'OCCIP', see
% classify_RAW_eeglab_data.m and classify_TFR_from_eeglab_data.m for more 
% options.
% gsettings.exclsubj = index numbers of subjects to skip (not including
% these subjects in the results structs). No subjects are excluded when
% left empty (default).
% Specify whether to extract the BDM or the FEM data using
% gsettings.plotmodel = 'BDM' (default) or 'FEM'
% You can also (optionally) specify which conditions you want to extract
% (and in which order) using gsettings.plotorder, as a cell array.
%
% Example: 
% gsettings.timelim = [0.2 1.2];
% gsettings.mpcompcor_method = 'cluster_based';
% gsettings.startdir = '/Volumes/backup/WM_debunk_EEG';
% gsettings.plotorder = { 'LOCATION_TASK','SEARCH_TASK'};
% gsettings.plotmodel = 'FEM';
% [stats, weights, gsettings] = compute_group_MVPA('',gsettings);
% 
% By J.J.Fahrenfort, VU 2015, 2016

% First get some settings
if nargin<3
    mask = [];
end

% Main routine, is a folder name specified? If not, pop up selection dialog
if isempty(folder_name)
    if ~isfield(gsettings,'startdir')
        gsettings.startdir = '';
        disp('NOTE: it is easier to select a directory when you indicate a starting directory in gsettings.startdir');
    end
    folder_name = uigetdir(gsettings.startdir,'select directory to plot');
    if ~ischar(folder_name)
        return
    end
    if ~isfield(gsettings,'plotorder') || isempty(gsettings.plotorder)
        dirz = dir(folder_name);
        dirz = {dirz([dirz(:).isdir]).name};
        dirz = dirz(cellfun(@isempty,strfind(dirz,'.')));
        % gsettings.plotorder = dirz;
    else
        dirz = gsettings.plotorder;
    end
    for cdirz = 1:numel(dirz)
        [stats(cdirz), weights(cdirz), gsettings] = subcompute_group_MVPA([folder_name filesep dirz{cdirz}],gsettings,mask);
    end
    gsettings.folder = folder_name;
else
    if ~exist('folder_name','dir') && ~iscell(folder_name) 
        error([folder_name ' should refer to a full and existing folder path']);
    end
    [stats, weights, gsettings] = subcompute_group_MVPA(folder_name,gsettings,mask);
end

function [stats,weights,gsettings] = subcompute_group_MVPA(folder_name,gsettings,mask)
% set defaults
one_two_tailed = 'two';
indiv_pval = .05;
cluster_pval = .05;
plotsubjects = false;
name = [];
channelpool = 'ALL';
mpcompcor_method = 'uncorrected';
plottype = '2D';
reduce_dims = [];
plotmodel = 'BDM';
timelim = [];
timelim2 = [];
trainlim = [];
testlim = [];
freqlim = [];
exclsubj= [];
get_dim = 'freq_time'; % default
% unpack graphsettings
v2struct(gsettings);
if strcmpi(plottype,'2D') && isempty(reduce_dims)
    reduce_dims = 'diag'; % 'avtrain' 'diag' 'avtest' 'avfreq'
    gsettings.reduce_dims = reduce_dims;
end
if strcmpi(get_dim,'frequency_time') || strcmpi(get_dim,'time_frequency') || strcmpi(get_dim,'time_freq')
    get_dim = 'freq_time';
    gsettings.get_dim = get_dim;
end

% general time limit
if ~isempty(timelim) % && isempty(trainlim) && isempty(testlim) timelim takes precedence
    trainlim = timelim;
    testlim = timelim;
    gsettings.trainlim = trainlim;
    gsettings.testlim = testlim;
end

% fill some empties
if isempty(channelpool)
    error('no channelpool was specified in settings, set gsettings.channelpool to ALL or OCCIP or similar, see classify_ functions.');
end
pval(1) = indiv_pval;
pval(2) = cluster_pval;

% some logical checking: is this a frequency folder?
f_time_time = ~isempty(dir([folder_name filesep channelpool filesep 'freq*'])); % contains_2Dfreqs
f_freq_time = ~isempty(dir([folder_name filesep channelpool filesep 'allfreqs'])); % contains_3Dfreqs
freqfolder = any([f_time_time f_freq_time]);
if freqfolder
    if strcmpi(get_dim,'freq_time') && f_freq_time
        plotFreq = {[filesep 'allfreqs']};
    elseif strcmpi(get_dim,'freq_time')
        error('ERROR: freq_time is not available in this folder, try gsettings.get_dim = time_time')
    else
        if isempty(freqlim)
            freqlim = input('What frequency or frequency range should I extract? ');
        end
        if numel(freqlim) == 1
            freqlim(2) = freqlim(1);
        end
    end
    if strcmpi(get_dim,'time_time') && f_time_time
        if numel(freqlim) > 1 % make a list of frequencies over which to average
            freqlist = dir([folder_name filesep channelpool filesep 'freq*']); freqlist = {freqlist(:).name};
            for c=1:numel(freqlist); freqsindir(c) = string2double(regexprep(freqlist{c},'freq','')); end
            freqsindir = sort(freqsindir); freqs2keep = find(freqsindir >= min(freqlim) & freqsindir <= max(freqlim));
            for c=1:numel(freqs2keep); plotFreq{c} = [filesep 'freq' num2str(freqsindir(freqs2keep(c)))]; end
        else % or just a single frequency
            plotFreq{1} = [filesep 'freq' num2str(freqlim)];
        end
    elseif strcmpi(get_dim,'time_time')
        disp('WARNING: time_time is not availaible in this folder, defaulting to freq_time')
        get_dim = 'freq_time';
        plotFreq{1} = [filesep 'allfreqs'];
    end
else
    plotFreq{1} = '';
    freqlim = [];
end
        
% if contains_2Dfreqs && ~isempty(freqlim)
%     if numel(freqlim) > 1 % make a list of frequencies over which to average
%         freqlist = dir([folder_name filesep channelpool filesep 'freq*']); freqlist = {freqlist(:).name};
%         for c=1:numel(freqlist); freqsindir(c) = string2double(regexprep(freqlist{c},'freq','')); end
%         freqsindir = sort(freqsindir); freqs2keep = find(freqsindir >= min(freqlim) & freqsindir <= max(freqlim));
%         for c=1:numel(freqs2keep); plotFreq{c} = [filesep 'freq' num2str(freqsindir(freqs2keep(c)))]; end
%     else % or just a single frequency
%         plotFreq{1} = [filesep 'freq' num2str(freqlim)];
%     end
% elseif contains_3Dfreqs
%     plotFreq{1} = [filesep 'allfreqs'];
% else
%     plotFreq{1} = '';
%     freqlim = [];
% end

% pack graphsettings with defaults
nameOfStruct2Update = 'gsettings';
gsettings = v2struct(freqlim,plotFreq,trainlim,testlim,one_two_tailed,indiv_pval,cluster_pval,plotmodel,mpcompcor_method,plottype,reduce_dims,freqlim,nameOfStruct2Update);

% make mask 2D if plot is 2D
if strcmpi(plottype,'2D')
    mask = diag(mask);
end

% get filenames
subjectfiles = dir([folder_name filesep channelpool plotFreq{1} filesep '*.mat']);
[~, condname] = fileparts(folder_name);
subjectfiles = { subjectfiles(:).name };

% limiting subjects
if ~isempty(exclsubj)
   subjectfiles = select_subjects(subjectfiles,exclsubj,true);
end

nSubj = numel(subjectfiles);
if nSubj == 0
    error(['cannot find data in specified folder ' folder_name filesep channelpool plotFreq]);
end
if plotsubjects;
    fh = figure('name','individual subjects');
    set(fh, 'Position', get(0,'Screensize'));
    set(fh,'color','w');
end

% do the loop, restrict time and frequency if applicable
for cSubj = 1:nSubj
    fprintf(1,'loading subject %d of %d\n', cSubj, nSubj);
    clear ClassOverTimeAv;
    % loop over frequencies
    for cFreq = 1:numel(plotFreq)
        matObj = matfile([folder_name filesep channelpool plotFreq{cFreq} filesep subjectfiles{cSubj}]);
        if ~isempty(whos(matObj,'BDM'))
            if strcmpi(plotmodel,'BDM')
                v2struct(matObj.BDM); % unpack fields
            else
                v2struct(matObj.FEM); % unpack fields
            end
        else % BW compatibility, will become obsolete over time
            ClassOverTime = matObj.ClassOverTime;
            WeightsOverTime = matObj.WeightsOverTime;
            if isfield(matObj,'covPatternsOverTime')
                covPatternsOverTime = matObj.covPatternsOverTime;
                corPatternsOverTime = matObj.corPatternsOverTime;
            end
        end
        % sum up frequency accuracies to compute average over frequencies
        if ~exist('ClassOverTimeAv','var');
            ClassOverTimeAv = zeros(size(ClassOverTime));
        end
        ClassOverTimeAv = ClassOverTimeAv + ClassOverTime;
    end
    % by default it computes the average over frequencies when specifying
    % time_time (gsettings.reduce_dims = 'avfreq' is actually superfluous)
    ClassOverTime = ClassOverTimeAv / numel(plotFreq);
    settings = matObj.settings;
    [ClassOverTime, settings, gsettings, lim1,lim2] = restrict_ClassOverTime(ClassOverTime,settings,gsettings);
    v2struct(gsettings);
    ClassTotal{1}(cSubj,1:size(ClassOverTime,1),1:size(ClassOverTime,2)) = ClassOverTime;
    
    % first get the relevant electrodes and put them in the same order
    % as in chanlocdata
    if ~isfield(settings,'chanlocs') % if no chanlocdata exist in settings
        if ~exist('chanlocdata','var')
            chanlocdata = readlocs('plotting_1005.sfp','importmode','native');
        end
        [~, chanindex, dataindex] = intersect({chanlocdata(:).labels},settings.channels,'stable');
        chanlocs = chanlocdata(chanindex); % put all in the same order as imported locations 
    else % otherwise just extract from settings
        chanlocs = settings.chanlocs{1};
        if ~exist('firstchanlocs','var')
            firstchanlocs = chanlocs;
        end
        [~, ~, dataindex] = intersect({firstchanlocs(:).labels},{chanlocs(:).labels},'stable');
        chanlocs = firstchanlocs;
    end
    if numel(chanlocs) < numel(settings.channels)
        error('could not find location info for all channels');
    end
    
    % STILL NEED TO IMPLEMENT FREQUENCY AVERAGING OVER WeightsOverTime,
    % covPatternsOverTime, corPatternsOverTime etc 
    
    % extract weights and patterns: time x channel or freq x time x channel
    % make sure that all subject data are in the same electrode order by using dataindex
    if strcmpi(settings.dimord,'frequency_time')
        if exist('WeightsOverTime', 'var'); WeightsOverTimeAll(cSubj,:,:,:,:) = WeightsOverTime(lim2,lim1,dataindex,:); end;
        if exist('covPatternsOverTime', 'var'); covPatternsOverTimeAll(cSubj,:,:,:) = covPatternsOverTime(lim2,lim1,dataindex); end;
        if exist('corPatternsOverTime', 'var'); corPatternsOverTimeAll(cSubj,:,:,:) = corPatternsOverTime(lim2,lim1,dataindex); end;
    else
        if exist('WeightsOverTime', 'var') && ~isempty(WeightsOverTime); WeightsOverTimeAll(cSubj,:,:,:) = WeightsOverTime(lim1,dataindex,:); end;
        if exist('covPatternsOverTime', 'var'); covPatternsOverTimeAll(cSubj,:,:) = covPatternsOverTime(lim1,dataindex); end;
        if exist('corPatternsOverTime', 'var'); corPatternsOverTimeAll(cSubj,:,:) = corPatternsOverTime(lim1,dataindex); end;
    end
    if exist('C2_average', 'var'); C2_averageAll(cSubj,:,:,:) = C2_average(lim2,lim1,:); end;
    if exist('C2_percondition', 'var'); C2_perconditionAll(cSubj,:,:,:,:) = C2_percondition(lim2,lim1,:,:); end;
    
    indivClassOverTime(cSubj,:,:) = ClassOverTime;
    if strcmpi(plottype,'3D')
        indivClassAv(cSubj) = mean(mean(ClassOverTime));
    else
        indivClassAv(cSubj) = mean(mean(diag(squeeze(ClassOverTime))));
    end
    if plotsubjects
        subplot(numSubplots(nSubj,1),numSubplots(nSubj,2),cSubj);
        onestat.ClassOverTime = ClassOverTime;
        onestat.StdError = [];
        onestat.pVals = [];
        onestat.indivClassOverTime = [];
        onestat.indivClassAv = [];
        onestat.settings = settings;
        onestat.condname = condname;
        onestat.channelpool = channelpool;
        tempsettings = gsettings;
        tempsettings.timetick = 500;
        plot_MVPA(onestat,tempsettings);
    end
end

% determine chance level
if strcmpi(settings.measuremethod,'hr-far') || strcmpi(plotmodel,'FEM')
    chance = 0;
else
    chance = 1/settings.nconds;
end

% compute standard errors and averages
ClassStdErr(1:size(ClassTotal{1},2),1:size(ClassTotal{1},3)) = std(ClassTotal{1},0,1)/sqrt(size(ClassTotal{1},1));
ClassAverage(1:size(ClassTotal{1},2),1:size(ClassTotal{1},3)) = mean(ClassTotal{1},1);
ClassTotal{2} = repmat(chance,size(ClassTotal{1}));

% statistical testing
if strcmpi(mpcompcor_method,'fdr')
    % FDR CORRECTION
    if strcmpi(one_two_tailed,'two')
        [~,ClassPvals(1:size(ClassTotal{1},2),1:size(ClassTotal{1},3))] = ttest(ClassTotal{1},ClassTotal{2},'tail','both');
    else
        [~,ClassPvals(1:size(ClassTotal{1},2),1:size(ClassTotal{1},3))] = ttest(ClassTotal{1},ClassTotal{2},'tail','right');
    end
    if strcmpi(plottype,'2D')
        thresh = fdr(diag(squeeze(ClassPvals)),pval(2));
    else
        thresh = fdr(squeeze(ClassPvals),pval(2));
    end
    ClassPvals(ClassPvals>thresh) = 1;
elseif strcmpi(mpcompcor_method,'cluster_based')
    % CLUSTER BASED CORRECTION
    if strcmpi(plottype,'2D') 
        for cSubj = 1:nSubj
            DiagTotal{1}(cSubj,:) = diag(squeeze(ClassTotal{1}(cSubj,:,:)));
            DiagTotal{2}(cSubj,:) = diag(squeeze(ClassTotal{2}(cSubj,:,:)));
        end
        % mask = diag(mask);
        [ DiagPvals, pStruct ] = cluster_based_permutation(DiagTotal{1},DiagTotal{2},gsettings,settings,mask);
        ClassPvals = ones(size(ClassAverage));
        ClassPvals(logical(eye(size(ClassPvals)))) = DiagPvals;
    else
        [ClassPvals(1:size(ClassTotal{1},2),1:size(ClassTotal{1},3)), pStruct] = cluster_based_permutation(ClassTotal{1},ClassTotal{2},gsettings,settings,mask);
    end
elseif strcmpi(mpcompcor_method,'uncorrected')
    % NO MP CORRECTION
    if strcmpi(one_two_tailed,'two')
        [~,ClassPvals(1:size(ClassTotal{1},2),1:size(ClassTotal{1},3))] = ttest(ClassTotal{1},ClassTotal{2},'tail','both');
    else
        [~,ClassPvals(1:size(ClassTotal{1},2),1:size(ClassTotal{1},3))] = ttest(ClassTotal{1},ClassTotal{2},'tail','right');
    end
    ClassPvals(~mask) = 1;
else
    % NO TESTING, PLOT ALL
    ClassPvals = zeros([size(ClassTotal{1},2) size(ClassTotal{1},3)]);
end

% outputs
stats.ClassOverTime = ClassAverage;
stats.StdError = ClassStdErr;
stats.pVals = ClassPvals;
stats.indivClassOverTime = indivClassOverTime;
stats.indivClassAv = indivClassAv;
stats.settings = settings;
stats.condname = condname;
stats.filenames = subjectfiles;
stats.channelpool = channelpool;
if exist('pStruct','var')
    stats.pStruct = pStruct;
end
%gsettings = v2struct(name,nameOfStruct2Update);

% compute weights stuff
weights.chanlocs = chanlocs;
if exist('WeightsOverTimeAll','var')
    weights.avWeights = squeeze(mean(WeightsOverTimeAll,1));
    weights.indivWeights = squeeze(WeightsOverTimeAll);
end
if exist('covPatternsOverTimeAll','var')
    weights.avCovPatterns = squeeze(mean(covPatternsOverTimeAll,1));   
    weights.indivCovPatterns = squeeze(covPatternsOverTimeAll);
end
if exist('corPatternsOverTimeAll','var')
    weights.avCorPatterns = squeeze(mean(corPatternsOverTimeAll,1));
    weights.indivCorPatterns = squeeze(corPatternsOverTimeAll);
end
if exist('C2_averageAll','var')
    weights.CTF = squeeze(mean(C2_averageAll,1));
    weights.semCTF = squeeze(std(C2_averageAll,0,1)/sqrt(size(C2_averageAll,1)));
    weights.indivCTF = squeeze(C2_averageAll);
end
if exist('C2_perconditionAll','var')
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
weights.filenames = subjectfiles;

% determine file name for graph when plotting
base_folder = folder_name;
if exist('startdir','var') && strncmp(base_folder,startdir,numel(startdir))
    tmp = base_folder(numel(startdir)+1:end);
    tmp = strsplit(tmp,filesep);
    try
        gsettings.outfile = strjoin('_',tmp(3:end-1));
    catch
        gsettings.outfile = strjoin(tmp(3:end-1),'_');
    end
else
    gsettings.outfile = '';
end
disp('done!');

% find limits within which to constrain ClassOverTime
function [limitedClassOverTime, settings, gsettings, lim1, lim2] = restrict_ClassOverTime(ClassOverTime,settings,gsettings) 
% unpack graphsettings
v2struct(gsettings);
v2struct(settings);
if strcmpi(dimord,'frequency_time')
    freqs = settings.freqs; % little hack to overcome the fact that freqs is also a function and v2struct can't deal with that
end

% OLD bw compatible, can remove self/other after some time
if isfield(settings.times,'self')
    clear times;
    times{1} = settings.times.self;
end
if isfield(settings.times,'other')
    times{2} = settings.times.other;
elseif numel(settings.times) == 1 && strcmpi(settings.dimord,'time_time')
    times{2} = times{1};
end

% continue limit operation
if strcmpi(dimord,'frequency_time') && numel(freqlim)>1
    lim2 = nearest(freqs,freqlim(1)):nearest(freqs,freqlim(2));
    settings.freqs = freqs(lim2);
elseif strcmpi(dimord,'frequency_time') && numel(freqlim) == 1
    lim2 = nearest(freqs,freqlim);
    settings.freqs = freqs(lim2);
elseif strcmpi(dimord,'time_time') && ~isempty(testlim)
    lim2 = nearest(times{2}*1000,testlim(1)):nearest(times{2}*1000,testlim(2));
    times{2} = times{2}(lim2);
else
    lim2 = true(size(ClassOverTime,1),1);
end

% the time dimension (always present)
if ~isempty(trainlim)
    lim1 = nearest(times{1}*1000,trainlim(1)):nearest(times{1}*1000,trainlim(2));
    times{1} = times{1}(lim1);
else
    lim1 = true(size(ClassOverTime,2),1);
end

% if plottype is 2D and the domain is time_time, restrict lim2 by the same thing as lim1
if strcmpi(plottype,'2D') && strcmpi(dimord,'time_time')
    lim2 = lim1;
end

% now restrict
limitedClassOverTime = ClassOverTime(lim2,lim1);
settings.times = times;

% and if applicable, average
if strcmpi(reduce_dims,'avfreq') && strcmpi(dimord,'freq_time')
    limitedClassOverTime = mean(limitedClassOverTime,1);
elseif strcmpi(reduce_dims,'avtrain') && strcmpi(dimord,'time_time')
    limitedClassOverTime = mean(limitedClassOverTime,1);
elseif strcmpi(reduce_dims,'avtest') && strcmpi(dimord,'time_time')
    limitedClassOverTime = mean(limitedClassOverTime,2);
elseif strcmpi(reduce_dims,'diag')
    if sum(size(squeeze(limitedClassOverTime))==1) % is one of the dimensions 1?
        limitedClassOverTime = diag(limitedClassOverTime); % put the values back on the diagonal for consistency
    end
end