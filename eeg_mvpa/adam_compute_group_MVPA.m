function [stats,cfg] = adam_compute_group_MVPA(cfg,folder_name,mask)
% ADAM_COMPUTE_GROUP_MVPA extracts group-level classification data from single-subject MVPA results
% that were given by ADAM_MVPA_FIRSTLEVEL, and if requested performs a statistical test at the group
% level. The result is a structure that contains the group average, together with the
% subject-specific classification data, classifier weights, and (if requested) statistical outcomes,
% and can be used for plotting (as an input for ADAM_PLOT_MVPA and ADAM_PLOT_MVPA_weights).
%
% Use as:
%   stats = adam_compute_group_MVPA(cfg)
%
% The cfg (configuration) input structure should at least specify the path to the directory where
% the results are located:
%
%       cfg.startdir         = string specifiying a directory where the results of 
%                              ADAM_MVPA_FIRSTLEVEL are located; cfg.startdir can be left empty, but
%                              this requires you to navigate from your Matlab root folder to the
%                              desired folder every time. Even if you specifiy a directory here,
%                              ADAM will nonetheless always generate a pop-up folder-selection
%                              window, giving you the option to switch to another folder, or to
%                              confirm by clicking 'Open'.
%
% The following options are optional:
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
%       cfg.plot_dim         = 'time_time' (default); or 'freq_time'
%       cfg.reduce_dims      = [] (default); this will take all dimensions of the level-1 result; to 
%                              customize, you can specify: 'diag' (only take the diagonal if the
%                              level-1 was time-by-time generalization), 'avtrain' (average over a
%                              time window used for training the classifier and plot every tested
%                              time point), 'avtest' (average over a time window used for testing
%                              the classifier and plot every trained time point) or 'avfreq'
%                              (average over a frequency band of interest).
%       cfg.trainlim         = [int int]; is the time limits over which to constrain the training 
%                              data, in ms.
%       cfg.testlim          = [int int]; is the time limits over which to constrain the testing 
%                              data, in ms.
%       cfg.timelim          = [int int]; constrains trainlim and testlim at once (takes precedence
%                              over trainlim and testlim)
%       cfg.freqlim          = [int int]; is the frequency limits over which to constrain the 
%                              frequency dimension. NOTE: if you have time-by-time classification
%                              for multiple frequencies, and you leave freqlim empty, ADAM will ask
%                              you in the Command window to specify a frequency to select or a
%                              frequency range to average over. This is because level-1 time-by-time
%                              for a frequency range is stored in separate folder for each
%                              frequency. If at the level-1 analysis a diagonal approach is
%                              specified, classification is saved in one frequency-by-time matrix,
%                              in which case leaving freqlim empty results in all frequencies
%                              present in the data being analyzed.
%       cfg.channelpool      = 'ALL_NOSELECTION' or other string, e.g. 'OCCIP', according to which
%                              channelpools have been analyzed at the level-1; only one pool can be
%                              specified per group-level analysis.
%       cfg.plot_model       = 'BDM' (default) or 'FEM'; Specify whether to extract the backward
%                              decoding (BDM) or forward encoding (FEM) data
%       cfg.plot_order       = string e.g. {'cond1','cond2'} to specify which conditions you want to 
%                              extract (and in which order).
%       cfg.plotsubjects     = false (default); or true; if true, during importing single-subject
%                              data, one figure with subplots is generated that displays each
%                              single-subject classification result.
%       cfg.exclsubj         = '' (default); index numbers of subjects to skip (not including these 
%                              subjects in the results structs). No subjects are excluded when left
%                              empty (default).
%
% The output stats structure will contain the following fields:
%
%       stats.ClassOverTime:        NxM matrix; group-average classification accuracy over N 
%                                   testing time points and M training time points; note that if
%                                   reduce_dims is specified, M will be 1, and ClassOverTime
%                                   will be squeezed to a Nx1 matrix of classification over time.
%       stats.indivClassOverTime:   PxNxM matrix; classificition accuracy over testing and training 
%                                   time for P subjects 
%       stats.StdError:             NxM matrix; standard-deviation across subjects over time-time
%       stats.pVals:                NxM matrix; p-values of each tested time-time point
%       stats.pStruct:              struct; cluster info, if mpcompcor_method was set to
%                                   'cluster_based'
%       stats.mpcompcor_method:     string; correction method ('uncorrected' is default)
%       stats.settings:             struct; the settings grabbed from the level-1 results
%       stats.condname:             string; name of the level-1 folder
%       stats.channelpool:          string; the selected channel pool
%       stats.weights:              struct; classification weights: group-average and
%                                   subject-specific, for actual weights and the
%                                   correlation/covariance class separability maps
%       stats.cfg:                  struct; the cfg used to create these stats
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Example usage: 
%
% cfg.startdir          = '/Volumes/project1/decoding/EEG/';
% cfg.timelim           = [-200 1000];
% cfg.mpcompcor_method  = 'cluster_based';
% cfg.plot_order        = { 'famous_nonfamous','famous_scrambled'};
% cfg.reduce_dims       = 'avtrain';
% cfg.trainlim          = [50 150];
%
% stats = adam_compute_group_MVPA(cfg);
% 
% part of the ADAM toolbox, by J.J.Fahrenfort, VU, 2017/2018
% 
% See also ADAM_COMPUTE_GROUP_ERP, ADAM_MVPA_FIRSTLEVEL, ADAM_PLOT_MVPA, ADAM_PLOT_BDM_WEIGHTS

% First get some settings
if nargin<3
    mask = [];
end
if nargin<2
    folder_name = '';
end

plot_order = {};
reduce_dims = '';
freqlim = [];

% backwards compatibility
plot_dim = 'time_time'; % default, 'time_time' or 'freq_time'
v2struct(cfg);
if exist('one_two_tailed','var')
    error('The cfg.one_two_tailed field has been replaced by the cfg.tail field. Please replace cfg.one_two_tailed with cfg.tail using ''both'', ''left'' or ''right''. See help for further info.');
end
if exist('plotmodel','var')
    plot_model = plotmodel;
    cfg.plot_model = plot_model;
    cfg = rmfield(cfg,'plotmodel');
end
if exist('get_dim','var')
    plot_dim = get_dim;
    cfg = rmfield(cfg,'get_dim');
end
if strcmpi(plot_dim,'frequency_time') || strcmpi(plot_dim,'time_frequency') || strcmpi(plot_dim,'time_freq')
    plot_dim = 'freq_time';
end
if exist('plotorder','var')
    plot_order = plotorder;
    cfg.plot_order = plot_order;
    cfg = rmfield(cfg,'plotorder');
end
cfg.plot_dim = plot_dim;

% check freqlimits
if (strcmpi(plot_dim,'freq_time') && ~strcmpi(reduce_dims,'avfreq')) && (numel(freqlim) == 1 || (~isempty(freqlim) && abs(diff(freqlim)) <= 2))
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
            [stats, cfg] = subcompute_group_MVPA(cfg, [folder_name filesep plot_order{cdirz}], mask);
        else % getting from multiple folders
            [stats(cdirz), cfg] = subcompute_group_MVPA(cfg, [folder_name filesep plot_order{cdirz}], mask);
        end
    end
else
    if ~iscell(folder_name)
        folder_name = {folder_name};
    end
    for cdirz = 1:numel(folder_name)
        if ~exist(folder_name{cdirz},'dir')
            error([folder_name ' should refer to a full and existing folder path. Alternatively leave folder_name empty to pop up a selection dialog.']);
        end
        [stats(cdirz), cfg] = subcompute_group_MVPA(cfg, [folder_name{cdirz}], mask);
    end
end

% subroutine for each condition
function [stats,cfg] = subcompute_group_MVPA(cfg, folder_name,mask)
% set defaults
tail = 'both';
indiv_pval = .05;
cluster_pval = .05;
plotsubjects = false;
name = [];
channelpool = ''; 
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
if isempty(channelpool)
    chandirz = dir(folder_name);
    chandirz = {chandirz([chandirz(:).isdir]).name};
    chandirz = sort(chandirz(cellfun(@isempty,strfind(chandirz,'.'))));
    channelpool = chandirz{1};
    disp(['No cfg.channelpool specified, defaulting to channelpool ' channelpool ]);
end

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
        % You can comment out the error line below if you are not fond of this type of error handling :-)
        error('You specified cfg.plot_dim = ''time_time'', but time_time is not available in this folder.');
        disp('WARNING: You specified cfg.plot_dim = ''time_time'', but time_time is not available in this folder. Defaulting to cfg.plot_dim = ''freq_time''');
        plot_dim = 'freq_time';
        plotFreq{1} = [filesep 'allfreqs'];
    end
else
    % You can comment out the error lines below if you are not fond of this type of error handling :-)
    if ~isempty(freqlim)
        error('You indicated a frequency or frequency range to extract using cfg.freqlim, but this is not a frequency results folder! Please select a folder that contains frequency results.');
    end
    if strcmpi(plot_dim,'freq_time')
         error('You indicated you wanted to plot frquency by time by specifying cfg.plot_dim = ''freq_time'', but this is not a frequency results folder! Please select a folder that contains frequency results.');
    end
    plotFreq = {''};
    freqlim = [];
end
        
% pack graphsettings with defaults
nameOfStruct2Update = 'cfg';
cfg = v2struct(freqlim,plotFreq,trainlim,testlim,tail,indiv_pval,cluster_pval,plot_model,mpcompcor_method,reduce_dims,freqlim,nameOfStruct2Update);

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
    error(['cannot find data in specified folder ' folder_name filesep channelpool plotFreq{:} ' maybe you should specify (a different) cfg.channelpool?']);
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
        onestat.cfg = [];
        tmpcfg = cfg;
        tmpcfg.plot_model = plot_model;
        tmpcfg.plotsubjects = true;
        tmpcfg.plot_order = {condname};
        adam_plot_MVPA(tmpcfg,onestat);
        subjname = subjectfiles{cSubj};
        underscores = strfind(subjname,'_');
        subjname = regexprep(subjname(underscores(2)+1:underscores(end)-1),'_',' ');
        ntitle(subjname,'fontsize',10,'fontweight','bold');
        drawnow;
    end
end

% determine chance level
if strcmpi(settings.measuremethod,'hr-far') || strcmpi(plot_model,'FEM')
    chance = 0;
else
    chance = 1/settings.nconds;
end

% compute standard errors and averages
%ClassStdErr(1:size(ClassOverTimeAll{1},2),1:size(ClassOverTimeAll{1},3)) = std(ClassOverTimeAll{1},0,1)/sqrt(size(ClassOverTimeAll{1},1));
ClassStdErr = shiftdim(std(ClassOverTimeAll{1},0,1)/sqrt(size(ClassOverTimeAll{1},1)));
if sum(sum(ClassStdErr)) == 0 ClassStdErr = []; end % don't plot stderror when there is none
%ClassAverage(1:size(ClassOverTimeAll{1},2),1:size(ClassOverTimeAll{1},3)) = mean(ClassOverTimeAll{1},1);
ClassAverage = shiftdim(mean(ClassOverTimeAll{1},1));
ClassOverTimeAll{2} = repmat(chance,size(ClassOverTimeAll{1}));

% statistical testing
if nSubj > 1
    if strcmpi(mpcompcor_method,'fdr')
        % FDR CORRECTION
        [~,ClassPvals(1:size(ClassOverTimeAll{1},2),1:size(ClassOverTimeAll{1},3))] = ttest(ClassOverTimeAll{1},ClassOverTimeAll{2},pval(1),tail);
        thresh = fdr(squeeze(ClassPvals),pval(2));
        ClassPvals(ClassPvals>thresh) = 1;
    elseif strcmpi(mpcompcor_method,'cluster_based')
        % CLUSTER BASED CORRECTION
        [ClassPvals, pStruct] = cluster_based_permutation(ClassOverTimeAll{1},ClassOverTimeAll{2},cfg,settings,mask);
    elseif strcmpi(mpcompcor_method,'uncorrected')
        % NO MP CORRECTION
        [~,ClassPvals(1:size(ClassOverTimeAll{1},2),1:size(ClassOverTimeAll{1},3))] = ttest(ClassOverTimeAll{1},ClassOverTimeAll{2},pval(1),tail);
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
if isfield(stats.cfg,'plotsubjects')
     stats.cfg = rmfield(stats.cfg,'plotsubjects');
end
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
if ~isfield(settings,'chanlocs') || isempty(settings.chanlocs{1}) % if no chanlocdata exist in settings
    if ~exist('chanlocdata','var')
        % chanlocdata = readlocs('plotting_1005.sfp','importmode','native');
        chanlocdata = readlocs(findcapfile,'importmode','native');
    end
    [~, chanindex, dataindex] = intersect({chanlocdata(:).labels},settings.channels{1},'stable');
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
        error('Cannot find data, select different location in the directory hierarchy and/or check path settings.');
    end
    folder_name = fullfile(folder_name,nextlevel{1});
    ndirs = ndirs + 1;
    containsmat = ~isempty(dir(fullfile(folder_name, '*.mat')));
    containsfreq = ~isempty(dir(fullfile(folder_name, 'freq*'))) || ~isempty(dir(fullfile(folder_name, 'allfreqs')));
    if containsmat || containsfreq
        notfound = false;
    end
end