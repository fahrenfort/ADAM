function stats = adam_compute_group_MVPA(cfg,folder_name)
% ADAM_COMPUTE_GROUP_MVPA computes group-level classification data from single-subject results that
% were computed by adam_MVPA_firstlevel. It also performs statistical tests at the group level.
% When executing adam_compute_group_MVPA, a selection dialog will pop up. This dialog allows the
% user to select a directory containing the single subject results (output of adam_MVPA_firstlevel)
% for which to compute the group stats variable. The user can either select a directory referring to
% a specific analysis (e.g. EEG_FAM_VS_NONFAMOUS), or select one directory higher up in the
% hierarchy (e.g. RAW_EEG) which may contains several such analyses. The function creates a stats
% output variable. This variable is a structure that contains the group average, together with the
% subject-specific classification data, classifier weights, and (if requested) statistical outcomes,
% and can be used as input for the plotting functions (e.g. adam_plot_MVPA and
% adam_plot_BDM_weights). When selecting a directory containing several analyses, the stats variable
% will be an array of which each element contains a group analysis.
%
% Use as:
%   stats = adam_compute_group_MVPA(cfg)
%
% The cfg (configuration) input structure can contain the following optional parameters:
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
%       cfg.iterations       = 1000 (default); the number of permutations used for cluster based
%                              permutation testing.
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
%       cfg.mask            =  Optionally, you can provide a mask: a binary matrix to pre-select a
%                              'region of interest' to constrain the comparison. You can for example
%                              base the mask on pVals matrix in stats.pVals (beware of double
%                              dipping).
%       cfg.plot_dim         = 'time_time' (default); or 'freq_time' when plotting time-frequency
%       cfg.reduce_dims      = '' (default); this will take all dimensions of the level-1 result; to 
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
%       cfg.symmetry         = false; (default). When set to true, the temporal generalization
%                              matrix is averaged across the diagonal, so that train to testing time
%                              generalizations and test to train time generalizations are
%                              effectively averaged. This can be a useful operation when train and
%                              test time differ in their signal to noise ratio, so that certain
%                              asymetries in the temporal generalization matrix are due to SNR
%                              differences rather than interesting asymetries in generalization. For
%                              example, some effects may show up when you train on time point X and
%                              test on time point Y, but not when you train on time point Y but not
%                              on time point X. Setting cfg.symmetry = true; resolves this by
%                              simply averaging X->Y and Y->X generalizations.
%       cfg.channelpool      = 'ALL_NOSELECTION' or other string, e.g. 'OCCIP', according to the
%                              electrode selection at the level-1; only one pool can be specified
%                              per group-level analysis. See the help of adam_MVPA_firstlevel and
%                              select_channels for details.
%       cfg.plot_model       = 'BDM' (default) or 'FEM'; Specify whether to extract the backward
%                              decoding (BDM) or forward encoding (FEM) results
%       cfg.plotsubjects     = false (default); or true; if true, the individual subject
%                              classification results of all subjects that are extracted will be
%                              plotted in a single figure, with a separate subplot for each subject
%                              to enable inspection of the data underlying group averages.
%       cfg.exclsubj         = {} (default); Cell array of strings containing the names of subjects
%                              to exclude from the group-level analysis, e.g. cfg.exclsubj =
%                              {'S01_faces_exp', 'S02_faces_exp'};. You can also use part of a name,
%                              e.g. it is sufficient to specify {'S01', 'S02'}; Subject names are
%                              based on the file names that were used as input for the first-level
%                              analysis.
%       cfg.read_confidence   = false (default) or true. If true, adam_compute_group_MVPA will
%                              attempt to read the confidence scores of the classifier for
%                              individual trials / time points / frequencies for every subject and
%                              store it in the stats output variable: stats.indivConf. The indivConf
%                              field will contain the subfields: indivConf(nSubj).scores,
%                              indivConf(nSubj).true_class_labels, indivConf(nSubj).event_labels,
%                              indivConf(nSubj).trial_index, indivConf(nSubj).dimord. When
%                              confidence scores are read into the stats variable, they can 
%                              be used by adam_correlate_confidence to produce a new
%                              stats variable that can subsequently be plotted. 
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
%       stats.latencies             struct; containing three fields with onset latencies, based on
%                                   the percent-amplitude latency (50%): single subject, jackknife
%                                   and GA (grand average). Note that single subject outcomes are
%                                   not very reliable, it is much better to rely on the jackknife
%                                   results. For details about latency computation see: Liesefeld,
%                                   H. R. (2018). Estimating the Timing of Cognitive Operations With
%                                   MEG/EEG Latency Measures: A Primer, a Brief Tutorial, and an
%                                   Implementation of Various Methods. Frontiers in Neuroscience,
%                                   12, 765. http://doi.org/10.3389/fnins.2018.00765
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
% cfg.startdir          = '/Volumes/project1/FACE_EXPERIMENT/RAW_EEG_RESULTS/';
% cfg.timelim           = [-200 1000];
% cfg.mpcompcor_method  = 'cluster_based';
% cfg.reduce_dims       = 'avtrain';
% cfg.trainlim          = [50 150];
% cfg.exclsubj          = {'S04', 'S07'};
%
% stats = adam_compute_group_MVPA(cfg);
% 
% part of the ADAM toolbox, by J.J.Fahrenfort, VU, 2017/2018
% 
% See also: ADAM_MVPA_FIRSTLEVEL, ADAM_COMPUTE_GROUP_ERP, ADAM_PLOT_MVPA, ADAM_PLOT_BDM_WEIGHTS, FDR_BH

% First get some settings
if nargin<2
    folder_name = '';
end
mask = [];
plot_order = {};
reduce_dims = '';
freqlim = [];
plotsubjects = false;
singleplot = false;
symmetry = false;

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
if exist('channels','var') && ~isfield(cfg,'channelpool')
    channelpool = channels;
    cfg.channelpool = channelpool;
end
cfg.singleplot = singleplot;

% check freqlimits
if (strcmpi(plot_dim,'freq_time') && ~strcmpi(reduce_dims,'avfreq')) && (numel(freqlim) == 1 || (~isempty(freqlim) && abs(diff(freqlim)) <= 2))
    disp('WARNING: your cfg.freqlim indicates a rather small range of frequencies given cfg.plot_dim ''freq_time'', use cfg.reduce_dims = ''avfreq'' if you intend to average. Now simply plotting all frequencies.');
    freqlim = [];
    cfg.freqlim = [];
end

% Main routine, is a folder name specified? If not, pop up selection dialog
if isempty(folder_name)
    if ~isfield(cfg,'startdir')
        cfg.startdir = '';
        disp('NOTE: it is easier to select a directory when you indicate a starting directory using cfg.startdir, otherwise you have to start selection from root every time...');
    end
    folder_name = uigetdir(cfg.startdir,'select directory for which to compute group MVPA results');
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

% when plotting individual subjects in a single plot also need to create a figure and set colors
if plotsubjects && singleplot
    fh = figure('name','individual subjects');
    set(fh, 'Position', get(0,'Screensize'));
    set(fh,'color','w');
    % set line colors
    if ~exist('line_colors','var')
        line_colors = {[.5 0 0] [0 .5 0] [0 0 .5] [.5 .5 0] [0 .5 .5] [.5 0 .5] [.75 0 0] [0 .75 0] [0 0 .75] [.75 .75 0] [0 .75 .75] [.75 0 .75] };
    end
    if ~iscell(line_colors)
        for c=1:size(line_colors,1)
            cell_colors{c} = line_colors(c,:);
        end
        line_colors = cell_colors;
    end
end

% loop through directories (results folders)
for cdirz = 1:numel(plot_order)
    if plotsubjects && singleplot
        cfg.line_color = line_colors{cdirz};
    end
    [stats(cdirz), cfg] = subcompute_group_MVPA(cfg, [folder_name filesep plot_order{cdirz}], mask);
end

% subroutine for each condition
function [stats,cfg] = subcompute_group_MVPA(cfg, folder_name,mask)
% set defaults
tail = 'both';
indiv_pval = .05;
cluster_pval = .05;
plotsubjects = false;
channelpool = ''; 
mpcompcor_method = 'uncorrected';
plot_model = 'BDM'; % 'BDM' or 'FEM'
reduce_dims = []; % 'diag' 'avtrain' 'avtest' or 'avfreq'
timelim = [];
trainlim = [];
testlim = [];
freqlim = [];
exclsubj = [];
compute_randperm = false;
read_confidence = false;
symmetry = false;
trysymmetry = false;
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
if isempty(channelpool)
    chandirz = dir(folder_name);
    chandirz = {chandirz([chandirz(:).isdir]).name};
    chandirz = sort(chandirz(cellfun(@isempty,strfind(chandirz,'.'))));
    channelpool = chandirz{1};
    disp(['No cfg.channelpool specified, defaulting to channelpool ' channelpool ]);
end

% some logical checking: is this a frequency folder?
freqfolder_contains_time_time = ~isempty(dir([folder_name filesep regexprep(channelpool,',','_') filesep 'freq*'])); 
freqfolder_contains_freq_time = ~isempty(dir([folder_name filesep regexprep(channelpool,',','_') filesep 'allfreqs']));
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
            freqlist = dir([folder_name filesep  regexprep(channelpool,',','_') filesep 'freq*']); freqlist = {freqlist(:).name};
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
        
% pack cfg with defaults
nameOfStruct2Update = 'cfg';
cfg = v2struct(freqlim,plotFreq,trainlim,testlim,tail,indiv_pval,cluster_pval,plot_model,mpcompcor_method,reduce_dims,freqlim,symmetry,nameOfStruct2Update);

% get filenames
subjectfiles = dir([folder_name filesep  regexprep(channelpool,',','_') plotFreq{1} filesep '*.mat']);
[~, condname, ext] = fileparts(folder_name);
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
    error(['cannot find data in specified folder ' folder_name filesep regexprep(channelpool,',','_') plotFreq{1} ' maybe you should specify (a different) cfg.channelpool?']);
end

% prepare figure in case individual subjects are plotted
if plotsubjects && ~singleplot
    fh = figure('name',['individual subjects: ' condname]);
    set(fh, 'Position', get(0,'Screensize'));
    set(fh,'color','w');
end

% do the loop, restrict time and frequency if applicable
trialcount = zeros(1,nSubj);
for cSubj = 1:nSubj
    fprintf(1,'loading subject %d of %d\n', cSubj, nSubj);
    
    % initialize subject
    clear ClassOverTimeAv WeightsOverTimeAv covPatternsOverTimeAv corPatternsOverTimeAv C2_averageAv C2_perconditionAv;

    % loop over frequencies (if no frequencies exist, it loads raw data)
    for cFreq = 1:numel(plotFreq)
        
        % locate data
        matObj = matfile([folder_name filesep regexprep(channelpool,',','_') plotFreq{cFreq} filesep subjectfiles{cSubj}]);
        if isempty(whos(matObj,'settings'))
            error([ subjectfiles{cSubj} ' in ' condname ' does not contain the expected fields, which suggests that the analysis did not complete succesfully. Check the error message that appeared during first level analysis.']);
        else
            settings = matObj.settings;
        end
        
        % for backward compatibility
        if strcmpi(settings.dimord,'frequency_time')
            settings.dimord = 'freq_time';
        end
        v2struct(settings);
        
        % compute trial counts
        trialcount(cSubj) = numel(trialinfo{1,1}) + numel(trialinfo{1,2});
        
        % get data
        if ~isempty(whos(matObj,'BDM')) && strcmpi(plot_model,'BDM')
            v2struct(matObj.BDM); % unpack
            if read_confidence
                if isempty(whos(matObj,'BDM_CONF')) % for backwards compatibility
                    BDM_CONF = 'No confidence scores of the classifier were saved during first level analysis. To save confidence scores, specify cfg.save_confidence = ''yes'' when running adam_MVPA_firstlevel.';
                else
                    BDM_CONF = matObj.BDM_CONF;
                end
                if ischar(BDM_CONF)
                    disp(BDM_CONF);
                    read_confidence = false; % nothing to read
                    cfg.read_confidence = read_confidence;
                end
            end
        elseif ~isempty(whos(matObj,'FEM')) && strcmpi(plot_model,'FEM')
            v2struct(matObj.FEM); % unpack
            % little hack to fix the measure method
            measuremethod = 'CTF slope';
            settings.measuremethod = measuremethod;
            read_confidence = false; % FEM does not have confidence scores (yet)
            cfg.read_confidence = read_confidence;
        else
            try % backward compatibility, can be removed in later version
                disp('WARNING: These seem to be analyses that were performed using ancient scripts. Unless this is pure replication, it is recommended to run the first levels again.');
                matObj = load([folder_name filesep regexprep(channelpool,',','_') plotFreq{cFreq} filesep subjectfiles{cSubj}]);
                v2struct(matObj);
                if ~exist('covPatternsOverTime', 'var')
                    try
                        covPatternsOverTime = PatternsOverTime;
                        corPatternsOverTime = PatternsOverTime;
                    catch
                        covPatternsOverTime = [];
                        corPatternsOverTime = [];
                    end
                end
                clear matObj PatternsOverTime;
            catch
                error('cannot find data');
            end
        end
        
        % find permutations and load them
        permfolder = [folder_name filesep regexprep(channelpool,',','_') plotFreq{cFreq} filesep 'randperm'];
        if compute_randperm && exist(permfolder,'dir')
            subjectpermutes = dir([permfolder filesep subjectfiles{cSubj}(1:end-4) '_PERM*.mat']);
            subjectpermutes = {subjectpermutes(:).name};
            countP = zeros(size(ClassOverTime));
            failed = 0;
            for cPerm=1:numel(subjectpermutes)
                if ~mod(cPerm,round(numel(subjectpermutes)/10))
                    fprintf(1,'reading random permutation %d of %d\n', cPerm, numel(subjectpermutes));
                end
                try
                    if strcmpi(plot_model,'BDM')
                        load([permfolder filesep subjectpermutes{cPerm}],'BDM');
                        countP = countP + (BDM.ClassOverTime >= ClassOverTime);
                    elseif strcmpi(plot_model,'FEM')
                        load([permfolder filesep subjectpermutes{cPerm}],'FEM');
                        countP = countP + (FEM.ClassOverTime >= ClassOverTime);
                    end
                catch ME
                    disp(ME.message);
                    failed = failed + 1;
                end
            end
            pVals = countP / (numel(subjectpermutes)-failed);
        end
        
        % initialize chanlocs
        if ~exist('firstchanlocs','var')
            firstchanlocs = [];
        end
        
        % find limits
        % [settings, cfg, lim1, lim2, dataindex, firstchanlocs, chanlocdata] = find_limits(settings, cfg, firstchanlocs, chanlocdata);
        [settings, cfg, lim1, lim2, dataindex ] = find_limits(settings, cfg, firstchanlocs);
        v2struct(cfg);
        
        % create symmetry across the diagonal?
        if symmetry
            if size(ClassOverTime,1) == size(ClassOverTime,2)
                tempClassOverTime = ClassOverTime;
                pairs = nchoosek(1:size(ClassOverTime,1),2);
                for cPairs = 1:size(pairs,1)
                    newPerf = mean([tempClassOverTime(pairs(cPairs,1),pairs(cPairs,2)) tempClassOverTime(pairs(cPairs,2),pairs(cPairs,1))]);
                    ClassOverTime(pairs(cPairs,1),pairs(cPairs,2)) = newPerf;
                    ClassOverTime(pairs(cPairs,2),pairs(cPairs,1)) = newPerf;
                end
                clear tempClassOverTime;
            elseif trysymmetry == false
                trysymmetry = true;
                disp('Cannot compute symmetrical temporal generalization as train and test axes are not the same size, will try again after limit operation.');
            end
        end
        
        % limit ClassOverTime
        ClassOverTime = ClassOverTime(lim1,lim2);
        
        % create symmetry across the diagonal?
        if symmetry && trysymmetry
            if size(ClassOverTime,1) == size(ClassOverTime,2)
                tempClassOverTime = ClassOverTime;
                pairs = nchoosek(1:size(ClassOverTime,1),2);
                for cPairs = 1:size(pairs,1)
                    newPerf = mean([tempClassOverTime(pairs(cPairs,1),pairs(cPairs,2)) tempClassOverTime(pairs(cPairs,2),pairs(cPairs,1))]);
                    ClassOverTime(pairs(cPairs,1),pairs(cPairs,2)) = newPerf;
                    ClassOverTime(pairs(cPairs,2),pairs(cPairs,1)) = newPerf;
                end
                clear tempClassOverTime;
            else
                disp('Cannot compute symmetrical temporal generalization as train and test axes are not the same size even after limit operation.');
            end
        end
        
        % limit Confidence
        if read_confidence
            if strcmpi(BDM_CONF.dimord,'rpt_freq_time_class')
                BDM_CONF.scores = BDM_CONF.scores(:,lim1,lim2,:);
            elseif strcmpi(BDM_CONF.dimord,'rpt_time_class')
                BDM_CONF.scores = BDM_CONF.scores(:,lim1,:);
            end
        end
        
        % hack to compute onset latencies for first level random permutations
        if compute_randperm && exist(permfolder,'dir')
            pVals = pVals(lim1,lim2);
            pVals = diag(pVals);
            
            % compute actual
            [obsPosSizes, obsNegSizes, posLabels, negLabels] = compute_clusters(diag(ClassOverTime),pVals,indiv_pval);

            % mask containing significant points
            sigMatrix = zeros(size(pVals));
            sigMatrix(pVals < indiv_pval) = 1;
            
            % iterate to find max clusters under random permutation
            % IMPLEMENT LATER WHEN I HAVE TIME
            
            % cut out insignificant points
            randSizes = 5; % quick hack to identify max clusters without iteration, simply cut off at 5
            labels = 1:max(unique(posLabels));
            sizes2rm = find(obsPosSizes < max(randSizes));
            for c = 1:numel(sizes2rm)
                sigMatrix(posLabels==labels(sizes2rm(c))) = 0;
                posLabels(posLabels==labels(sizes2rm(c))) = 0;
            end
            pVals = ~sigMatrix;
            pStruct.posclusters = compute_pstruct(posLabels,pVals,diag(ClassOverTime),cfg,settings);
%             obsolete now:
%             if ~isempty(pStruct.posclusters)
%                 onsetTimes(cSubj) = min([ pStruct.posclusters(:).start_time ]);
%             else
%                 onsetTimes(cSubj) = nan;
%             end
        end
        
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
            if isempty(testlim)
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
        if read_confidence
            indivConf(cSubj) = BDM_CONF;
        end
    else
        indx = [{cSubj} repmat({':'}, 1, ndims(C2_average))];
        C2_averageAll(indx{:}) = C2_average;
        indx = [{cSubj} repmat({':'}, 1, ndims(C2_percondition))];
        C2_perconditionAll(indx{:}) = C2_percondition;
    end
    
    % plot individual subjects
    if plotsubjects
        if singleplot
            hold on; % for plotting multiple lines
        end
        subplot(numSubplots(nSubj,1),numSubplots(nSubj,2),cSubj);
        onestat.ClassOverTime = ClassOverTime;
        onestat.StdError = [];
        if exist('pVals','var')
            onestat.pVals = pVals;
        else
            onestat.pVals = [];
        end
        onestat.indivClassOverTime = [];
        onestat.settings = settings;
        onestat.condname = condname;
        onestat.channelpool = channelpool;
        onestat.cfg = [];
        tmpcfg = cfg;
        if singleplot
            tmpcfg.line_colors = {line_color}; % for plotting multiple lines
        end
        tmpcfg.plot_model = plot_model;
        tmpcfg.plotsubjects = true;
        tmpcfg.plotsigline_method = 'both';
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
if any(strcmpi(settings.measuremethod,{'hr-far','dprime','hr','far','mr','cr'})) || strcmpi(settings.measuremethod,'\muV') || ~isempty(strfind(settings.measuremethod,' difference')) || strcmpi(plot_model,'FEM')
    chance = 0;
elseif strcmpi(settings.measuremethod,'AUC')
    chance = .5;
else
    chance = 1/settings.nconds;
end
settings.chance = chance;

% compute standard errors and averages
ClassStdErr = shiftdim(squeeze(std(ClassOverTimeAll{1},0,1)/sqrt(size(ClassOverTimeAll{1},1))));
if sum(sum(ClassStdErr)) == 0 ClassStdErr = []; end % don't plot stderror when there is none
ClassAverage = shiftdim(squeeze(mean(ClassOverTimeAll{1},1)));
ClassOverTimeAll{2} = repmat(chance,size(ClassOverTimeAll{1}));

% statistical testing
% note that 'alpha' and 'tail' do not need to be specified explictly, the ttest and ttest2 functions
% are backwards compatible with Matlab 2012b
if nSubj > 1
    if strcmpi(mpcompcor_method,'fdr')
        % FDR CORRECTION
        [~,ClassPvals] = ttest(ClassOverTimeAll{1},ClassOverTimeAll{2},indiv_pval,tail);
        ClassPvals = squeeze(ClassPvals);
        pValsUncorrected = ClassPvals;
        h = fdr_bh(ClassPvals,cluster_pval,'dep');
        ClassPvals(~h) = 1;
        pStruct = compute_pstructs(h,ClassPvals,ClassOverTimeAll{1},ClassOverTimeAll{2},cfg,settings);
    elseif strcmpi(mpcompcor_method,'cluster_based')
        % CLUSTER BASED CORRECTION
        [~,ClassPvals] = ttest(ClassOverTimeAll{1},ClassOverTimeAll{2},indiv_pval,tail);
        ClassPvals = squeeze(ClassPvals);
        pValsUncorrected = ClassPvals;
        [ClassPvals, pStruct] = cluster_based_permutation(ClassOverTimeAll{1},ClassOverTimeAll{2},cfg,settings,mask);
    elseif strcmpi(mpcompcor_method,'uncorrected')
        % NO MP CORRECTION
        [h,ClassPvals] = ttest(ClassOverTimeAll{1},ClassOverTimeAll{2},indiv_pval,tail);
        ClassPvals = squeeze(ClassPvals);
        pStruct = compute_pstructs(squeeze(h),ClassPvals,ClassOverTimeAll{1},ClassOverTimeAll{2},cfg,settings);
    else
        % NO TESTING, PLOT ALL
        ClassPvals = [];
    end
else
    if compute_randperm
        ClassPvals = pVals;
    else
        ClassPvals = [];
    end
end
ClassPvals = shiftdim(squeeze(ClassPvals));

% outputs
stats.ClassOverTime      = ClassAverage;
stats.StdError           = ClassStdErr;
stats.pVals              = ClassPvals;
if exist('pValsUncorrected','var')
    stats.pValsUncorrected   = pValsUncorrected;
end
stats.mpcompcor_method   = mpcompcor_method;
stats.indivClassOverTime = ClassOverTimeAll{1};
if read_confidence
    stats.indivConf      = indivConf;
end
stats.settings           = settings;
stats.trialcount         = trialcount;
stats.condname           = condname;
stats.filenames          = subjectfiles;
stats.channelpool        = channelpool;
if exist('pStruct','var')
    stats.pStruct        = pStruct;
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
% compute latency
try
    stats.latencies = extract_latency(cfg,stats);
catch ME
    disp('Cannot extract latencies.');
    disp(ME.message);
    stats.latencies = [];
end
% obsolete now:
% if exist('onsetTimes')
%     stats.onsetTimes = onsetTimes;
% end

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

function [settings, cfg, lim1, lim2, dataindex] = find_limits(settings, cfg, firstchanlocs) 
% find limits within which to constrain ClassOverTime
times = []; % need to initialize this, because times also happens to be a variable
v2struct(cfg); % unpack cfg
v2struct(settings); % unpack settings
if strcmpi(dimord,'freq_time') 
    freqs = settings.freqs; % to fix that freqs is also a function and v2struct can't deal with that
end

% some backwards compatibility
if numel(times) == 1 && strcmpi(settings.dimord,'time_time')
    if iscell(times)
        times{2} = times{1};
    elseif isfield(times,'self')
        clear times;
        times{1} = settings.times.self;
        try
            times{2} = settings.times.other;
        catch
            times{2} = times{1};
        end
    end
end

% get the relevant electrodes and obtain the correct order for weights (from the first subject)
% extract channels from settings
channels = settings.channels;
if iscell(channels)
    if iscell(channels{1})
        channels = channels{1}; % take the training channel list
    end
end
% extract the channel locations from settings
if isfield(settings,'chanlocs')
    chanlocs = settings.chanlocs;
    if iscell(chanlocs)
        chanlocs = chanlocs{1};
    end
else
    % if no chanlocdata exist, read them in
    if exist('1005chanlocdata.mat','file')
        load('1005chanlocdata.mat');
    else
        chanlocdata = readlocs(trycapfile,'importmode','native'); % from standard 10-20 system
    end
    [~, ~, chanindex] = intersect(channels,{chanlocdata(:).labels},'stable'); % takes FT_EEG.label as point of departure!
    chanlocs = chanlocdata(chanindex);
    if numel(chanlocs) ~= numel(channels)
        error('Could not read in all channel positions for some reason. Either make sure they exist in the original EEG file, or add these positions to the 1005chanlocdata.mat file so they can be obtained during processing.');
    end 
end

% make sure that channel locations and channel labels are in the same order
if sum(strcmpi({chanlocs(:).labels}, channels')) ~= numel(channels)
    error('The order of the channel labels and the channel locations is different for some reason, revisit your original EEG files to see what might be wrong.');
end 
% put channels and channel locations in the same order as the first channel locations that were read in
if isempty(firstchanlocs)
    firstchanlocs = chanlocs;
end
[~, ~, dataindex] = intersect({firstchanlocs(:).labels},{chanlocs(:).labels},'stable');
settings.channels = channels(dataindex);
settings.chanlocs = chanlocs(dataindex);

% when plotting diagonal, time axes need to be the same between train and test
if strcmpi(reduce_dims,'diag') && strcmpi(dimord,'time_time')
     if ~isempty(trainlim) && ~isempty(testlim)
         trainlim(1) = max([trainlim(1) testlim(1)]);
         trainlim(2) = min([trainlim(2) testlim(2)]);
         testlim = trainlim;
     else
         if ~isempty(trainlim)
             testlim = trainlim;
         elseif ~isempty(testlim)
             trainlim = testlim;
         else
             % taking the maximal overlapping time window between train and test
             trainlim(1) = max([min(times{1}) min(times{2})])*1000;
             trainlim(2) = min([max(times{1}) max(times{2})])*1000;
             testlim = trainlim;
         end
     end
end

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