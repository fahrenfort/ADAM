function classify_TFR_from_eeglab_data(filepath,filenames,outpath,nFolds,channelset,method,crossclass_and_or_resample,tf_and_erp_baseline,frequencies,varargin)
% classify_TFR_from_eeglab_data is an internal function of the ADAM toolbox to do TFR extraction and
% MVPA classification in one step. Refer to the help of ADAM_MVPA_FIRSTLEVEL for proper instructions
% on how to use the ADAM toolbox.
% classify_TFR_from_eeglab_data uses compute_TFR_from_eeglab.m, see help for this  function to learn
% more about the TFR input parameters to method (which can also contain values like 'total',
% 'induced' or 'evoked'). filenames either contains a single filename for testing and training, or
% two filenames separated by a comma (the first for training, second for testing) outpath contains
% the folder where results should be stored (if empty defaults to filepath) nFolds is the number of
% cohorts in which the data is split up for the leave-one-out procedure (only when the same data are
% used for training and testing). If there are separate sets for testing and training, nFolds
% defaults to 1. In that case, testing and training need to contain the same event codes, or you
% need to specify separate event codes for the testing set by separating by a semicolon (see below).
% channelset determines the electrode subset that is used for testing. This can be done either
% numerically (assumes a 64 channel 10-20 system):
% 1 ALL (uses all electrodes)
% 2 OCCIP (uses occipital electrodes)
% 3 PARIET (uses parietal electrodes)
% 4 FRONTAL (uses frontal electrodes)
% 5 TEMPORAL (uses temporal electrodes)
% 6 OCCIPARIET (uses occipitoparietal electrodes)
% method specifies the method used for classification. Default is 'linear' for linear discriminant
% analysis. Other options are: 'diagLinear', 'mahalanobis' or 'quadratic'. 
% The dependent classification measure can be specified in method, either using 'AUC' (default),
% 'accuracy', or various SDT measures ('dprime', 'hr-far', etc). If more than two categories are
% present when computing an SDT measure, the script defaults back to computing AUC. Specify by
% adding to method, e.g. like this: method = 'linear,accuracy'.
% If you want to output the confidence scores that were assigned by the classifier to each of the
% trials rather than computing classifier performance, specify cfg.save_confidence = 'yes'.
% In method you can also specify that the classifier should be trained under random permutation of
% the training labels by adding 'randperm', separated by a comma, like this: method =
% 'linear,randperm'. In method you can also specify that the training labels should be permuted
% under random permutation by adding 'randperm', separated by a comma, like this: method =
% 'linear,randperm'. The results of computing under random permutation are stored in a folder called
% 'randperm' and can subsequently be used for hypothesis testing. Similarly, you can also run
% another complete iteration of the same data (all folds) without randomly permuting the labels, by
% specificying method = 'iterate'. The results of computing another iteration are stored in a folder
% called 'iterations' and can subsequently be used for averaging to get a cleaner result. You can
% specify the number of times you want to execute a permutation or iteration through the
% create_qsub_files function, by setting the number of repeats using settings.repeat as an argument
% to that function. In method, you can also specify whether the trial labels that are specified for
% each condition should be binned (averaged) to generate new trials, or not. This is done by adding
% 'bin' to the method, separated by comma's, like this: method = 'linear,bin'; If left out, no
% binning is applied (default). One can also specify 'bintrain' or 'bintest' if binning should only
% take place on the training or on the testing side. Other options for method are to specify whether
% to compute the Scalp Current Density (SCD) prior to classification. This is done by either
% specifying 'csd' or 'scd' (e.g. method = 'linear,csd') using the finite method from
% ft_scalpcurrentdensity from the fieldtrip toolbox. You can specify the number of times you want to
% execute a permutation or iteration through the create_qsub_files function, by setting the number
% of repeats using settings.repeat as an argument to that function. crossclass_and_or_resample
% specifies whether a cross-classification is performed over time, in which each time point is used
% to classify all other time points. This is very time consuming, better set to 0 for first
% analysis, but set to 1 if you want to get a better feel of the data. You can speed up cross
% classification by downsampling the data. This is done by specifying the new sampling rate,
% separating the value for whether or not to cross classify from the new sampling rate using a
% semicolon, like this: crossclass_and_or_resample = '1;128'; meaning that cross classification will
% be performed and that all data will be resampled to 128 Hz. If no value is specified to resample
% to, the data is left at the original sampling rate. frequencies allows you to select for which
% frequencies you want to compute classification accuracy. This can be specified as an array, as a
% list of comma separated values, or as a single value. If left empty or set to 0, all frequencies
% are classified. varargin is a variable set of conditions to select from, by which each main
% condition is specified as string with comma separated values:
% cond1 = '1,2,3';
% cond2 = '4,5,6';
% In that setup, all trials containing either a event code 1, 2 or 3 will be trained and tested as
% category 1, while all trials containing a 4, 5 and 6 will be trained and tested as category 2. If
% different condition nrs are used in training versus testing, separate the condition specificiation
% by a semicolon (can only be used when using different sets for testing and training):
% cond1 = '1,2;8,9';
% cond2 = '3;10';
% here 1, 2 (the first) would be for training and 8,9 (the second) for testing category 1, while 3
% would be training category 2 and 10 for testing category 2.
%
% J.J.Fahrenfort, VU 2014, 2016, 2018
% 
% Internal function of the ADAM toolbox by J.J.Fahrenfort, VU 2014, 2015, 2018
%
% See also: ADAM_MVPA_FIRSTLEVEL

% sanity checking and parameter extraction
warning('off','all')
if isempty(outpath)
    outpath = filepath;
else
    if ~exist(outpath,'dir')
        mkdir(outpath);
    end
end
if ~iscell(filenames)
    filenames = regexp(filenames, ';|,', 'split');
end
if numel(filenames)>2
    error('too many input files');
end
if isempty(nFolds)
    nFolds = 1;
end
if ischar(nFolds)
    nFolds = string2double(nFolds);
end
if isempty(method)
    method = 'linear';
end
if ischar(frequencies)
    frequencies = str2num(frequencies);
end
if isempty(crossclass_and_or_resample)
    crossclass = false;
    resample_eeg = false;
end
if ischar(crossclass_and_or_resample)
    crossclass_and_or_resample = string2double(crossclass_and_or_resample);
end
if numel(crossclass_and_or_resample) > 1
    crossclass = crossclass_and_or_resample(1);
    resample_eeg = crossclass_and_or_resample(2);
else
    crossclass = crossclass_and_or_resample;
    resample_eeg = 0;
end
if ischar(tf_and_erp_baseline) && ~isempty(tf_and_erp_baseline)
    tf_and_erp_baseline = string2double(tf_and_erp_baseline);
    if isempty(tf_and_erp_baseline)
        error('you did not specify the baseline(s) correctly, something went wrong during conversion');
    end
end
if size(tf_and_erp_baseline,1) == 2
   tf_baseline = tf_and_erp_baseline(1,:);
   erp_baseline = tf_and_erp_baseline(2,:);
else
   tf_baseline = tf_and_erp_baseline;
   erp_baseline = tf_and_erp_baseline;
end
if isempty(tf_baseline) || all(tf_baseline==0) || any(isnan(tf_baseline))
    tf_baseline = 'no';
end
if isempty(erp_baseline) || all(erp_baseline==0) || any(isnan(erp_baseline))
    erp_baseline = 'no';
end

% determine condition train-test classes
condSet = [];
if exist('varargin','var') && numel(varargin) > 0 && iscell(varargin)
    if ischar(varargin{1})
        for cCond = 1:numel(varargin)
            condSet{cCond} = string2double(varargin{cCond});
        end
    elseif isnumeric(varargin{1})
        condSet = varargin;
    end
end

% determine settings/methods for analysis
orig_method = method;
setmethod = method;
methods = regexp(method,',','split');
method = 'linear';
measuremethod = 'AUC';
randomize_labels = false;
iterate = false;
bintest = false;
bintrain = false;
compute_performance = true;
save_confidence = false;
tfr_method = 'total';
do_csd = false;
use_splines = false;
test_total = false;
do_FEM = false;
do_BDM = false;
clean_muscle = false;
clean_window = [];
basis_sigma = 1; % default width of basis set if not a delta, if this is empty, do a simple basis set (delta function)
unbalance_events = false;
unbalance_classes = false;
whiten = false;
whiten_test_using_train = false;
reproduce = false;
for c=1:numel(methods)
    if any(strcmpi(methods{c},{'linear', 'quadratic', 'diagLinear', 'diagQuadratic', 'mahalanobis'})) == 1
        method = methods{c};
    end
    if any(strcmpi(methods{c},{'total', 'evoked', 'induced'}))
        tfr_method = methods{c};
    end
    if strcmpi(methods{c},'test_total')
        test_total = true;
    end
    if strcmpi(methods{c},'bin')
        bintrain = true;
        bintest = true;
    end
    if strcmpi(methods{c},'bintrain')
        bintrain = true;
    end
    if strcmpi(methods{c},'bintest')
        bintest = true;
    end
    if strcmpi(methods{c},'randperm')
        iterate = true;
        randomize_labels = true;
    end
    if strcmpi(methods{c},'iterate')
        iterate = true;
    end
    if any(strcmpi(methods{c},{'hr-far','dprime','hr','far','mr','cr','accuracy','AUC'}))
        measuremethod = methods{c};
        if numel(condSet) ~= 2 && any(strcmpi(measuremethod,{'hr-far','dprime','hr','far','mr','cr'}))
            disp('Number of stimulus classes is unequal to 2, defaulting to computing AUC');
            measuremethod = 'AUC'; % defaulting back to AUC
        end
    end
    if any(strcmpi(methods{c},{'no_performance'}))
        compute_performance = false;
    end
    if any(strcmpi(methods{c},{'save_confidence'}))
        save_confidence = true;
    end
    if any(strcmpi(methods{c},{'csd','scd'}))
        do_csd = true;
    end
    if any(strcmpi(methods{c},{'splines','spline'}))
        use_splines = true;
    end
    if any(strcmpi(methods{c},{'FEM','do_FEM','forward'}))
        do_FEM = true;
    end
    if any(strcmpi(methods{c},{'BDM','do_BDM','backward'}))
        do_BDM = true;
    end
    if strcmpi(methods{c},'FEM_simple')
        basis_sigma = 0; % the same as setting basis_sigma to 0
    end
    if strncmpi(methods{c},'clean',5)
        clean_muscle = true;
        if numel(methods{c}) > 5
            clean_window = string2double(methods{c}(6:end))';
            if isnan(clean_window)
                clean_window = [];
            end
        end
    end
    if strncmpi(methods{c},'sigma',5)
        if numel(methods{c}) > 5
            basis_sigma = string2double(methods{c}(6:end))';
            if isnan(basis_sigma)
                basis_sigma = 1; % defaulting back to 1
            end
        end
    end
    if any(strcmpi(methods{c},{'unbalance_triggers', 'unbalance_events','unbalance','unbalanced'}))
        unbalance_events = true;
    end
    if any(strcmpi(methods{c},{'unbalance_classes'}))
        unbalance_classes = true;
    end
    if any(strcmpi(methods{c},{'undersample'}))
        disp('WARNING: between-class undersampling has become obsolete, always oversampling using ADASYN');
        unbalance_classes = false;
    end
    if any(strcmpi(methods{c},{'oversample'}))
        unbalance_classes = false;
    end
    if any(strcmpi(methods{c},{'whiten'}))
        whiten = true;
    end
    if any(strcmpi(methods{c},{'nowhiten'}))
        whiten = false;
    end
    if any(strcmpi(methods{c},{'reproduce'}))
        reproduce = true;
    end
end
if ~do_FEM && ~do_BDM
    do_BDM = true;
end
if save_confidence && ~do_BDM
    do_BDM = true;
    disp('Saving confidence scores currently only works with BDM (update for FEM is pending)');
end
% check condset
if isempty(condSet)
    error('Cannot find usable event specification.');
end
% double condset for test and train if only one is specified
if size(condSet{1},1) == 1
    condSet = put_this_condset(condSet,condSet,2);
end
% are train and test class specifications overlapping between each other?
trainevents = cellfun(@unique,get_this_condset(condSet,1),'UniformOutput',false);
testevents = cellfun(@unique,get_this_condset(condSet,2),'UniformOutput',false);
overlapping = false;
allthesame = false;
if any(ismember(unique([trainevents{:}]),unique([testevents{:}])))
    overlapping = true;
end
if all(ismember(unique([trainevents{:}]),unique([testevents{:}])))
    allthesame = true;
end
% are class specifications overlapping within train or within test?
if numel(unique([trainevents{:}])) < numel([trainevents{:}]) || numel(unique([testevents{:}])) < numel([testevents{:}])
    error('ERROR: There is an overlap between events in your class definitions. Two different classes cannot contain the same event code.');
end
% check nFolds
if numel(filenames) > 1 && nFolds > 1
    wraptext('WARNING: You specified different filenames for training and testing, with more than 1 fold. Leave-one-out is not applicable here. Defaulting nFolds to 1.',80);
    nFolds = 1; % if you do want to cut up independent sets into folds, you can do so by turning this safety check off
end
% check if condsets are non-overlapping while nFolds > 1, if so lower nFolds to 1
if numel(filenames) == 1 && nFolds > 1 && ~overlapping
    wraptext('WARNING: You specified non-overlapping event codes for classes in training and testing, with more than 1 fold. Leave-one-out is not applicable here. Defaulting nFolds to 1.',80);
    nFolds = 1; % if you do want to cut up independent sets into folds, you can do so by turning this safety check off
end
% if using same events for training and testing, increase nFolds
if numel(filenames) == 1 && nFolds == 1 && overlapping
    nFolds = 10;
    wraptext('WARNING: You dirty double dipper! You are using the same data for testing and training without a leave-one-out procedure. Defaulting nFolds to 10 for crossvalidation.',80);
end
% check if condsets are not the same but overlapping while events are balanced, if so unbalance_events
if numel(filenames) == 1 && ~unbalance_events && ~allthesame && overlapping
    unbalance_events = true;
    wraptext('WARNING: Some stimulus events overlap between train and test, within class balancing has now been turned OFF');
end
% check whether compute_performance is used properly
if ~compute_performance && nFolds > 1
    compute_performance = true;
    wraptext('WARNING: Always computing performance when using leave-on-out. Only set cfg.compute_performance = ''no'' when you do not have useful class labels to compute performance with. Defaulting back to cfg.compute_performance = ''yes''.',80);
end
% display class specification
wraptext('These are the class specifications. Each row contains the event codes that go into a single class (first row training, second row testing):',80);
celldisp(condSet,'class_spec');
% set random number generator
if reproduce
    rng('default');
else
    rng('shuffle','twister');
end

% display class specification
wraptext('These are the class specifications. Each row contains the event codes that go into a single class (first row training, second row testing):',80);
celldisp(condSet,'class_spec');

% Determine bundle name and/or electrode selection
[channelset, bundlename_or_bundlelabels] = return_channel_bundle(channelset);

% create a folder for this electrode group
outpath = [outpath filesep channelset];
if ~exist(outpath,'dir')
    mkdir(outpath);
end

% load data 
for cFile = 1:numel(filenames)
    msettings = [];
    msettings.channelpool = bundlename_or_bundlelabels;
    msettings.erp_baseline = erp_baseline; % NOTE: different in RAW.
    msettings.resample_eeg = false; % NOTE: this line is different in RAW. No resampling done here (yet)...
    msettings.do_csd = do_csd;
    msettings.clean_data = clean_muscle;
    msettings.clean_window = clean_window;
    msettings.shuffle_trials = true;
    [FT_EEG(cFile), filenames{cFile}, chanlocs{cFile}] = read_raw_data(filepath,filenames{cFile},outpath,msettings);
    % randomize labels for first level random permutation testing. NOTE: permuting all observations/labels regardless of the conditions in the experiment
    if randomize_labels
        FT_EEG(cFile).trialinfo = FT_EEG(cFile).trialinfo(randperm(numel(FT_EEG(cFile).trialinfo)));
    end
end

% duplicate data for testing if only one file is available
if numel(filenames) == 1
    FT_EEG(2) = FT_EEG;
    chanlocs{2} = chanlocs{1};
end

% extract some relevant trial info, training and testing data
for cSet = 1:2
    thisCondSet = get_this_condset(condSet,cSet);
    if unbalance_events
        trialinfo{cSet} = FT_EEG(cSet).trialinfo;
    else
        % bin/balance dataset (default action, this is not to achieve actual binnning, it just applies within-class balancing of conditions)
        FT_EEG_BINNED(cSet) = compute_bins_on_FT_EEG(FT_EEG(cSet),thisCondSet,'trial','original');
        trialinfo{cSet} = FT_EEG_BINNED(cSet).trialinfo;
        oldindex{cSet} = FT_EEG_BINNED(cSet).oldindex;
        % a bit of hack to assign event number -99 to unbalanced events (i.e. by selecting those events from thisCondSet that were not used to compute FT_EEG_BINNED)
        FT_EEG(cSet).trialinfo(setdiff(find(ismember(FT_EEG(cSet).trialinfo,[thisCondSet{:}])),[oldindex{cSet}{:}])) = -99;
    end
    % compute ERPs (baseline corrected, resampled, channels already selected)
    FT_ERP{cSet} = compute_erp_on_FT_EEG(FT_EEG(cSet),thisCondSet,'trial','bin');
    % also compute TFR for entire set
    [~, FT_TFR{cSet}] = compute_TFR_from_FT_EEG(FT_EEG(cSet),thisCondSet,resample_eeg,orig_method,tf_baseline,erp_baseline,frequencies);
    % keep track of channels and time line
    channels{cSet} = FT_TFR{cSet}.label;
    times{cSet} = FT_TFR{cSet}.time;
end

% if testing and training are on different files, check consistency
if numel(filenames) > 1 && ~all(strcmpi(channels{1},channels{2}))
    error('The electrodes do not occur in the same order in testing and training, some coding required to fix this...');
end

% Generate indices for folds to do training and testing
[setindex{1}, setindex{2}, nFolds] = make_folds(trialinfo{1},trialinfo{2},condSet,nFolds,~compute_performance);

% create file name based on train-test procedure
if numel(filenames) == 1
    if nFolds == 1
        train = ['_train_' cell2csv(cellfun(@(x) cell2csv(x,false,'-'),get_this_condset(condSet,1),'UniformOutput',false),false,'_')];
        test = ['_test_' cell2csv(cellfun(@(x) cell2csv(x,false,'-'),get_this_condset(condSet,2),'UniformOutput',false),false,'_')];
        filename = ['CLASS_PERF_' filenames{1} train test ];
    else
        filename = ['CLASS_PERF_' filenames{1} '_' num2str(nFolds) 'fold'];
    end
else
    filename = ['CLASS_PERF_train_' filenames{1} '_test_' filenames{2}];
end
if numel(filename) > 255
    filename = filename(1:255);
end

% save ERPs and TFRs
if ~crossclass && ~iterate % do not save for cross classification or iterations, otherwise every frequency would contain and FT_ERP and FT_TFR
    % a folder for time by frequency
    fullpath = fullfile(outpath, 'allfreqs');
    mkdir(fullpath);
    fullfilename = fullfile(fullpath, filename);
    save(fullfilename, 'FT_ERP', 'FT_TFR', '-v7.3');
end
clear FT_EEG_BINNED FT_ERP FT_TFR; % save memory by clearing

% between-class balancing: balance class instances by oversampling (only applied to training set)
if unbalance_classes
    wraptext('Between-class balancing is OFF. Make sure you know what you are doing, this can have undesirable effects on classifier bias when you have unevevenly represented stimulus classes in your design.',80);
else
    wraptext('Between-class balancing is ON, so that each stimulus class is evenly represented in the training set. If this is undesirable behavior, specify ''unbalance_classes'' in your methods.',80);
end

% within-class balancing: balance events within classes
if unbalance_events
    wraptext('Within-class balancing is OFF, such that an unequal distribution of event codes is allowed to contribute to each stimulus class. Make sure you know what you are doing, this can have undesirable effects on how you interpret your results.',80);
else
    [setindex{1}, setindex{2}] = unpack_binned(setindex{1}, setindex{2}, oldindex{1}, oldindex{2}); % unpack setindex{1} and setindex{2} to get back the original index
    wraptext('Within-class balancing is ON, such that event codes are evenly represented within each stimulus class. If event codes are very unevenly represented in your data, this can result in the loss of many trials. It does however, enforce a balanced design, which is important for interpretation. If this is undesirable behavior, specify ''unbalance_events'' in your methods.',80);
end

% generate folds, compute some stuff on them and save them temporarily before running MVPA
% do trial selection prior to TFR computation (important for induced!)
settrialinfo = [];
settrialindex = [];
for cFld = 1:nFolds
    for cSet = 1:2
        % FYI
        set_tfr_method{cSet} = tfr_method;
        % select trials belonging to this subset
        trialindex = vertcat(setindex{cSet}{cFld,:});
        % keep trying if running into memory issues (often temporary because all analyses are running in parallel)
        success = false; counterr = 0;
        while ~success
            try
                FT_EEG_2use = select_trials_from_FT_EEG(FT_EEG(cSet),trialindex);
                success = true;
            catch ME
                disp(ME); disp('Memory problem? Let''s wait a bit before trying again.'); counterr = counterr + 1; if counterr > 10; error('errrrr, tried 10 times, giving up now...'); end; pause(600); 
            end
        end
        method_2use = orig_method;
        % get some goodies
        trialinfo{cSet} = FT_EEG_2use.trialinfo;
        origtrialindex{cSet} = FT_EEG(cSet).origindex(trialindex)';
        % FYI
        set_use_splines(cFld,cSet) = use_splines;
        % if test_total is true, refrain from using induced during test, otherwise check trial
        % numbers, if too few subtract splines rather than average ERP (default)
        if strcmpi(tfr_method,'induced') && cSet == 2 && test_total
            disp('Using total power instead of induced power for testing, useful when there are too few trials in test to compute quality ERPs.');
            method_2use = strrep(orig_method,'induced','total');
            set_tfr_method{cSet} = 'total'; % FYI
        elseif strcmpi(tfr_method,'induced') && numel(trialindex) < (25*numel(condSet)) % default action if there are not enough trials in test
            disp('WARNING: fewer than 25 trials to compute an ERP for subtraction in induced power, computing spline on ERP before subtraction to be safe.');
            method_2use = [orig_method ',splines'];
            set_use_splines(cFld,cSet) = true; % FYI
        end
        % get relevant condSet for TFR computation
        condSet_2use = get_this_condset(condSet,cSet);
        % compute time frequency power spectrum and save result for each fold
        % keep trying if running into memory issues (often temporary)
        success = false; counterr = 0;
        while ~success
            try
                TFR_foldsave = compute_TFR_from_FT_EEG(FT_EEG_2use,condSet_2use,resample_eeg,method_2use,tf_baseline,erp_baseline,frequencies);
                success = true;
            catch ME
                disp(ME); disp('Memory problem? Let''s wait a bit before trying again.'); counterr = counterr + 1; if counterr > 10; error('errrrr, tried 10 times, giving up now...'); end; pause(600); 
            end
        end
        % temporarily save folded data
        [~,tmpf,~] = fileparts(tempname); % generates random filename
        fnames{cFld,cSet} = [filepath filesep '..' filesep TFR_foldsave.fname '_' filenames{1} '_fold' num2str(cFld) '_' num2str(cSet) '_' tmpf  '.mat'];
        save(fnames{cFld,cSet},'-v7.3','-struct','TFR_foldsave');
        clear FT_EEG_2use TFR_foldsave;
        % this is a check to make sure that the file has written to disk before continueing, which seemed to cause problems
        while ~exist(fnames{cFld,cSet},'file'); end; filesize = 0; sizenow = 1000;
        while filesize ~= sizenow 
            info = dir(fnames{cFld,cSet});
            filesize = info.bytes;
            disp(['file is ' num2str(filesize) ' bytes...']);
            pause(.5);
            info = dir(fnames{cFld,cSet});
            sizenow = info.bytes;
        end
        % just FYI, how big are the temporary files
        filesizes_MB(cFld,cSet) = round(sizenow/(2^20)*100)/100;
    end
    % FYI, store for every fold
    settrialinfo = [settrialinfo; trialinfo];
    settrialindex = [settrialindex; origtrialindex];
end % end folds loop in which temp files are created

clear FT_EEG; % clear the dataset so we don't need it in memory during analyses

% if no frequencies are specified, or if this is not a crossclassification, compute all frequencies
if isempty(frequencies) || sum(frequencies) == 0 || ~crossclass
    matObj1 = matfile(fnames{1,1});
    frequencies = matObj1.freq;
end
frequencies = round(frequencies*100)/100;

% now loop over frequencies
for cFreq = 1:numel(frequencies)
          
    % and run models, with as little overhead as possible
    clear BDM_* FEM_*;

    for cFld=1:nFolds
        % create file pointers for this fold
        frequency = frequencies(cFreq);
        [matObj1, dim_params1] = read_mat_file(fnames{cFld,1},channels{1},frequency);
        [matObj2, dim_params2] = read_mat_file(fnames{cFld,2},channels{2},frequency);
               
        % FYI
        fprintf(1,['fold: ' num2str(cFld) ', frequency: ' num2str(frequency) '\n']);
        
        % load data and put in correct format for BDM_and_FEM_FT_EEG, put channels in correct order to be sure
        train_FT_EEG.trial = squeeze(matObj1.powspctrm(dim_params1.index{:}));
        train_FT_EEG.dimord = regexprep(matObj1.dimord,{'freq_', '_freq'},''); % keep original dimord but cut out frequency
        train_FT_EEG.trialinfo = matObj1.trialinfo;
        train_condSet = get_this_condset(condSet,1);
        
        % also for testing
        test_FT_EEG.trial = squeeze(matObj2.powspctrm(dim_params2.index{:}));
        test_FT_EEG.dimord = regexprep(matObj2.dimord,{'freq_', '_freq'},'');
        test_FT_EEG.trialinfo = matObj2.trialinfo;
        test_condSet = get_this_condset(condSet,2);
        
        % fix dimord for whiten_FT_EEG and balance_FT_EEG functions to be sure (seem correct already)
        train_FT_EEG = fix_dimord(train_FT_EEG,'rpt_chan_time'); % should be trial * channel * time
        test_FT_EEG = fix_dimord(test_FT_EEG,'rpt_chan_time'); % should be trial * channel * time
        
        % (1) whiten train and test data
        if whiten
            [train_FT_EEG, FT_IE] = whiten_FT_EEG(train_FT_EEG,train_condSet);
            if nFolds > 4 || whiten_test_using_train
                % use train covariance to pre-whiten test data
                whiten_test_using_train = true; % to keep track in settings
                [test_FT_EEG] = whiten_FT_EEG(test_FT_EEG,test_condSet,FT_IE);
                clear FT_IE;
            else
                % or re-compute covariance matrix if test data is fully independent or at least 25% of the total data
                clear FT_IE;
                [test_FT_EEG] = whiten_FT_EEG(test_FT_EEG,test_condSet,[]);
            end
        end
        % (2) between-class balance train by oversampling minority class using ADASYN/SMOTE
        % STILL NEED TO CHECK INSIDE balance_FT_EEG IF FT_EEG IS NOT ALREADY BINNED:
        if ~unbalance_classes && ~unbalance_events % only works when events within classes are balanced, otherwise just discard, not that important
            train_FT_EEG = balance_FT_EEG(train_FT_EEG,train_condSet,whiten);
        end
        
        % fix dimord for BDM_and_FEM_FT_EEG function, could be standardized to save memory peaks, but no biggie
        train_FT_EEG = fix_dimord(train_FT_EEG,'chan_time_rpt'); % should be channel * time * trial
        test_FT_EEG = fix_dimord(test_FT_EEG,'chan_time_rpt'); % should be channel * time * trial

        % settings for backward and/or forward modelling
        msettings = [];
        msettings.crossclass = crossclass;
        msettings.method = method;
        msettings.measuremethod= measuremethod;
        msettings.compute_performance = compute_performance;
        msettings.save_confidence = save_confidence;
        msettings.doBDM = do_BDM;
        msettings.doFEM = do_FEM;
        msettings.basis_sigma = basis_sigma;
        msettings.fname = fnames{cFld,1}; % trainfile to read in raw data for covpatterns
        msettings.dim_params = dim_params1; % plus relevant info
        
        % run BDM and FEM
        [BDM, FEM] = BDM_and_FEM_FT_EEG(train_FT_EEG,test_FT_EEG,train_condSet,test_condSet,msettings);
        
        % clear obsolete data 
        clear train_FT_EEG test_FT_EEG train_condSet test_condSet
        
        % some BDM and FEM specific stuff
        if do_BDM
            % keep track of Classification over time for every fold (this is NaN when compute_performance == false)
            BDM_ClassOverT(cFld,:,:) = BDM.ClassOverTime;
            % keep track of confidence scores of the classifier for every fold
            if save_confidence
                BDM_ScoresOverT{cFld} = BDM.ScoresOverTime;
                BDM_TrueClassLabels{cFld} = BDM.TrueClassLabels;
            end
            % keep track of weights etc for every fold
            BDM_WeightsOverT(cFld,:,:) = BDM.WeightsOverTime; % fld x time x elec
            BDM_covPatternsOverT(cFld,:,:) = BDM.covPatternsOverTime; % fld x time x elec
            BDM_corPatternsOverT(cFld,:,:) = BDM.corPatternsOverTime; % fld x time x elec
        end
        if do_FEM
            FEM_ClassOverT(cFld,:,:) = tuning_from_matrix(FEM.C2_average,'slope',crossclass); % fld x t1 x t2
            FEM_WeightsOverT(cFld,:,:,:) = FEM.WeightsOverTime; % fld x time x elec x channel_response
            FEM_C2_average(cFld,:,:,:) = FEM.C2_average; % fld x time x time x channel_response
            FEM_C2_percondition(cFld,:,:,:,:) = FEM.C2_percondition; % fld x time x time x cond x channel_response
        end
        
    end % end folds loop, only frequency loop is left
    
    % mean data, average over folds where applicable
    BDM = [];
    FEM = [];
    if do_BDM
        BDM.ClassOverTime = squeeze(mean(BDM_ClassOverT,1)); % t1 x t2
        BDM.WeightsOverTime = squeeze(mean(BDM_WeightsOverT,1)); % time x elec x channel_response
        BDM.covPatternsOverTime = squeeze(mean(BDM_covPatternsOverT,1)); % time x elec
        BDM.corPatternsOverTime = squeeze(mean(BDM_corPatternsOverT,1)); % time x elec
        % save confidence scores for all test trials, requires a little reshuffling
        if save_confidence
            orig_index = vertcat(settrialindex{:,2});   % contains the original indices of the test trials
            event_labels = vertcat(settrialinfo{:,2});  % contains the original labels of the test trials
            scores = vertcat(BDM_ScoresOverT{:});       % put them into a single array for all folds
            true_class_labels = vertcat(BDM_TrueClassLabels{:});
            [trial_index, sort_index ] = sort(orig_index);   % sort into the right order
            scores = scores(sort_index,:,:,:);
            event_labels = event_labels(sort_index);
            true_class_labels = true_class_labels(sort_index);
            BDM_CONF.scores = scores;
            BDM_CONF.true_class_labels = true_class_labels;
            BDM_CONF.event_labels = event_labels;
            BDM_CONF.trial_index = trial_index;
            BDM_CONF.dimord = 'rpt_time_class'; % this may later be extended to also allow traintime x testtime, but this requires lots of memory
        else
            BDM_CONF = 'No confidence scores of the classifier were saved during first level analysis. To save confidence scores, specify cfg.save_confidence = ''yes'' when running adam_MVPA_firstlevel.';
        end
        clearvars BDM_* -except BDM_CONF;
    end
    if do_FEM
        FEM.ClassOverTime = squeeze(mean(FEM_ClassOverT,1));
        FEM.WeightsOverTime = squeeze(mean(FEM_WeightsOverT,1));
        FEM.C2_average = squeeze(mean(FEM_C2_average,1));
        FEM.C2_percondition = squeeze(mean(FEM_C2_percondition,1));
        clear FEM_*;
        BDM_CONF = [];
    end
    
    % save some settings so we know wtf just happened
    settings.nconds = numel(condSet);
    settings.nfolds = nFolds;
    settings.filenames = filenames;
    settings.crossclass = crossclass;
    settings.erp_baseline = erp_baseline;
    settings.tf_baseline = tf_baseline;
    settings.clean_window = clean_window;
    settings.BDM = do_BDM;
    settings.FEM = do_FEM;
    settings.basis_set_sigma = basis_sigma;
    settings.method_string = setmethod;
    settings.tfr_method = set_tfr_method;
    settings.use_splines_for_erps = set_use_splines; % when induced this means it subtracted the splined erps, when evoked it this means it used the splines of the erps
    settings.channelset = channelset;
    settings.channels = channels;
    settings.chanlocs = chanlocs;
    settings.times = times;
    settings.measuremethod = measuremethod;
    settings.compute_performance = compute_performance;
    settings.save_confidence = save_confidence;
    settings.trialinfo = settrialinfo;
    settings.trialindex = settrialindex;
    settings.condset = condSet;
    settings.csd_transform = do_csd;
    settings.bintrain = bintrain;
    settings.bintest = bintest;
    settings.unbalance_events = unbalance_events;
    settings.unbalance_classes = unbalance_classes;
    settings.filesizes_MB = filesizes_MB;
    settings.whiten = whiten;
    settings.whiten_test_using_train = whiten_test_using_train;
    
    % SAVE RESULTS
    % if crossclass == true, save crossclassification result PER FREQUENCY
    if crossclass
        settings.frequency = frequency;
        settings.dimord = 'time_time';
        % a folder for this frequency
        freqpath = fullfile(outpath, ['freq' num2str(frequency)]);
        % count filenames from 0001 onwards if computing under random permutation or iteration
        % create a folder for iterations / random permutations
        if iterate && randomize_labels
            freqpath = fullfile(freqpath, 'randperm');
        elseif iterate
            freqpath = fullfile(freqpath, 'iterations');
        end
        % create folder if it does not exist
        if ~exist(freqpath,'dir')
            mkdir(freqpath);
        end
        % count filenames from 001 onwards if computing under permutation or iteration
        if iterate
            fullfilename = find_filename(freqpath,filename);
        else
            fullfilename = fullfile(freqpath,filename);
        end
        save(fullfilename, 'FEM', 'BDM', 'BDM_CONF', 'settings', '-v7.3'); % not saving ERP and TFR to save memory
    % TEMP STORE freq x time RESULTS
    % if crossclass == false, store diagonals in freq x time matrix and save outside freq loop
    elseif ~crossclass
        if do_BDM
            BDMfreq_ClassOverTime(cFreq,:) = diag(BDM.ClassOverTime); % result is freq x time
            BDMfreq_WeightsOverTime(cFreq,:,:) = BDM.WeightsOverTime; % result is freq x time x chan
            BDMfreq_covPatternsOverTime(cFreq,:,:) = BDM.covPatternsOverTime; % result is freq x time x chan
            BDMfreq_corPatternsOverTime(cFreq,:,:) = BDM.corPatternsOverTime; % result is freq x time x chan     
            % save confidence scores for every frequency
            if save_confidence
                BDMfreq_scores(:,cFreq,:,:) = BDM_CONF.scores;
                BDMfreq_dimord = 'rpt_freq_time_class';
            else
                BDMfreq_scores = 'To save confidence scores, specify cfg.save_confidence = ''yes''';
            end
        end
        if do_FEM
            FEMfreq_ClassOverTime(cFreq,:) = diag(FEM.ClassOverTime);
            FEMfreq_WeightsOverTime(cFreq,:,:,:) = FEM.WeightsOverTime;
            for cTime=1:size(FEMfreq_ClassOverTime,2) % slightly more complex, can't use diag
                FEMfreq_C2_percondition(cFreq,cTime,:,:) = squeeze(FEM.C2_percondition(cTime,cTime,:,:));
                FEMfreq_C2_average(cFreq,cTime,:) = squeeze(FEM.C2_average(cTime,cTime,:));
            end
        end
    end % endif crossclass conditional
    
end % end frequency loop

% all frequencies done, delete obsolete files
for cFld = 1:nFolds
    for cSet = 1:2
        delete(fnames{cFld,cSet});
    end
end

% if crossclass == false, SAVE freq x time RESULTS (outside freq loop)
if ~crossclass
    BDM = [];
    FEM = [];
    if do_BDM
        BDM.ClassOverTime = BDMfreq_ClassOverTime;
        BDM.WeightsOverTime = BDMfreq_WeightsOverTime;
        BDM.covPatternsOverTime = BDMfreq_covPatternsOverTime;
        BDM.corPatternsOverTime = BDMfreq_corPatternsOverTime;
        if save_confidence
            BDM_CONF.scores = BDMfreq_scores;
            BDM_CONF.dimord = BDMfreq_dimord;
        end
        clearvars BDMfreq_*;
    end
    if do_FEM
        FEM.ClassOverTime = FEMfreq_ClassOverTime;
        FEM.WeightsOverTime = FEMfreq_WeightsOverTime;
        FEM.C2_percondition = FEMfreq_C2_percondition;
        FEM.C2_average = FEMfreq_C2_average;
        clear FEMfreq_*;
        BDM_CONF = [];
    end
    settings.freqs = frequencies;
    settings.dimord = 'freq_time';
    % a folder for time by frequency
    fullpath = fullfile(outpath, 'allfreqs');
    % count filenames from 0001 onwards if computing under random permutation or iteration
    % create a folder for iterations / random permutations
    if iterate && randomize_labels
        fullpath = fullfile(fullpath, 'randperm');
    elseif iterate
        fullpath = fullfile(fullpath, 'iterations');
    end
    % create folder if it does not exist
    if ~exist(fullpath,'dir')
        mkdir(fullpath);
    end
    % determine filename and save
    if iterate
        fullfilename = find_filename(fullpath,filename);
        save(fullfilename, 'FEM', 'BDM', 'BDM_CONF', 'settings', '-v7.3'); % not saving ERP and TFR to save memory
    else
        fullfilename = fullfile(fullpath,filename);
        save(fullfilename, 'FEM', 'BDM', 'BDM_CONF', 'settings', '-v7.3', '-append'); % this file also contains the ERPs and TFRs, so append
    end
        
end

% turn warnings back on
warning('on','all');

function findfile = find_filename(path,filename)
c = 1;
findfile = fullfile(path, sprintf([filename '_PERM%04d'], c));
while numel(dir([findfile '.*']))>0
    c = c + 1;
    findfile = fullfile(path, sprintf([filename '_PERM%04d'], c));
end