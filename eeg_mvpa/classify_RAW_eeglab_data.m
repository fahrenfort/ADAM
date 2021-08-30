function classify_RAW_eeglab_data(filepath,filenames,outpath,nFolds,channelset,method,crossclass_and_or_resample,erp_baselines,varargin)
% classify_RAW_eeglab_data is an internal function of the ADAM toolbox. Loads eeglab or FT data (no
% time frequency decomposition) and performs a multivariate classification procedure.
% Refer to the help of ADAM_MVPA_FIRSTLEVEL for proper instructions on how to use the ADAM toolbox.
% Filenames either contains a single filename for testing and training, or two filenames separated
% by a semi-colon or comma (the first for training, second for testing). outpath contains the folder
% where results should be stored (if empty defaults to filepath). nFolds is the number of cohorts in
% which the data is split up for the leave-one-out procedure (only when the same data are used for
% training and testing of course). When nFolds > 1, the algorithm divides trials over the folds such
% that each fold uses a fraction of (nFolds-1)/nFolds of all trials for training, and uses a
% fraction of 1/nFolds of all trials for testing. The training and testing procedure is applied to
% each fold separately, and the results of each of the folds is averaged at the end, ensuring that
% all data is tested once, without ever using the same data for training and testing. If some - but
% not all - trial labels are the same in the training and testing set (e.g. when some labels are
% used for training and testing, but other labels are unique to the training and/or testing set),
% the same procedure of splitting up trials for testing and training is followed separately for
% trials that are common to training and testing and the ones that are unique to training and
% testing, ensuring that their relative contribution during testing and training remains the same as
% it is in the original dataset. If you are using the same  dataset for training and testing, only
% set nFolds to 1 if training and testing labels are completely independent in the dataset. If there
% are separate sets for testing and training, nFolds defaults to 1. channelset determines the
% electrode subset that is used for testing. This can be done either numerically (assumes a 64
% channel 10-20 system):
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
% If you want to output the confidence
% scores that were assigned by the classifier to each of the trials rather than computing classifier
% performance, specify cfg.save_confidence = 'yes';. 
% In method you can also specify that the classifier should be trained under random permutation of
% the training labels by adding 'randperm', separated by a comma, like this: method =
% 'linear,randperm'. the results of computing under random permutation are stored in a folder called
% 'randperm' and counted so that they can be used for non-parametric hypothesis testing. Similarly,
% you can also run another complete iteration of the same data (all folds) without randomly
% permuting the labels, by specificying method = 'iterate'. The results of computing another
% iteration are stored in a folder called 'iterations' and counted so that they can subsequently be
% used for averaging to get a cleaner result. You can specify the number of times you want to
% execute a permutation or iteration through the create_qsub_files function, by setting the number
% of repeats using settings.repeat as an argument to that function. In method, you can also specify
% whether the trial labels that are specified for each condition should be binned (averaged) to
% generate new trials, or not. This is done by adding 'bin' to the method, separated by comma's,
% like this: method = 'linear,bin';
% If left out, no binning is applied (default). One can also specify 'bintrain' or 'bintest' if
% binning should only take place on the training or on the testing side.
% Other options for method are to specify whether to compute the Scalp Current Density (SCD) prior
% to classification. This is done by either specifying 'csd' or 'scd' (e.g. method = 'linear,csd')
% using the finite method from ft_scalpcurrentdensity from the fieldtrip toolbox.
% crossclass_and_or_resample specifies whether a cross-classification is performed over time, in
% which each time point is used to classify all other time points. This is very time consuming,
% better set to 0 for first analysis, but set to 1 if you want to get a better feel of the data. You
% can speed up cross classification by downsampling the data. This is done by specifying the new
% sampling rate, separating the value for whether or not to cross classify from the new sampling
% rate using a semicolon, like this:
% crossclass_and_or_resample = '1;128'; meaning that cross classification will be performed and that
% all data will be resampled to 128 Hz. If no value is specified to resample to, the data is left at
% the original sampling rate.
% erp_baselines allows you to perform baseline correction on the data. Specify in seconds
% [begin,end]. If specified as 0 or [0,0], no baseline correction will be applied (default). It is
% also possible to specify separate baselines for two input sets when separating them using a
% semicolon, like this: erp_baselines = '-.25,0; 0,0', which would apply a baseline correction of
% -250 to 0 ms only to the training set only. varargin is a variable set containing the event
% codes to select conditions from, specified as string with comma separated values:
% cond1 = '1,2,3';
% cond2 = '4,5,6';
% In that setup, all trials containing either an event code 1, 2 or 3 will be taken as belonging to
% condition 1, while all trials containing a 4, 5 and 6 will be taken as belonging to condition 2.
% If different event codes are used in training versus testing, separate the condition
% specification by a semicolon:
% cond1 = '1,2;8,9';
% cond2 = '3;10';
% here 1, 2 (the first) would be for training and 8,9 (the second for testing) condition 1, while 3
% would be training condition 2 and 10 for testing condition 2.
%
% Usage examples:
% classify_RAW_data('/home/eeg/TFR','subject1data1,subject1data2','home/eeg/mvpa/TFR',1,1,'linear,bintrials,hr-far',0,'-.2,0;-.4,-.2','1,2,3','4,5,6')
% classify_RAW_data('/home/eeg/TFR','subject1data','home/eeg/mvpa/TFR',4,2,'linear',0,0,'1;3','2;4')
%
% Internal function of the ADAM toolbox by J.J.Fahrenfort, VU 2014, 2015, 2018
%
% See also: ADAM_MVPA_FIRSTLEVEL

% first some sanity checks and initializations
warning('off','all')
if nargin < 10
    error('this function requires at least 10 arguments')
end
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
    resample_eeg = false;
end
if ischar(erp_baselines)
    erp_baselines = string2double(erp_baselines);
end
if size(erp_baselines,1) == 2
    erp_baseline{1} = erp_baselines(1,:);
    erp_baseline{2} = erp_baselines(2,:);
else
    erp_baseline{1} = erp_baselines;
    erp_baseline{2} = erp_baselines;
end
if isempty(erp_baseline{1}) || all(erp_baseline{1}==0) || any(isnan(erp_baseline{1}))
    erp_baseline{1} = 'no';
end
if isempty(erp_baseline{2}) || all(erp_baseline{2}==0) || any(isnan(erp_baseline{2}))
    erp_baseline{2} = 'no';
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
setmethod = method;
methods = regexp(method,',','split');
method = 'linear';
subtr_method = 'subtr_bin';
measuremethod = 'AUC';
randomize_labels = false;
iterate = false;
bintest = false;
bintrain = false;
compute_performance = true;
save_confidence = false;
do_csd = false;
compute_induced = false;
do_FEM = false;
do_BDM = false;
clean_muscle = false;
clean_window = [];
basis_sigma = 1; % default width of basis set if not a delta, if this is empty, do a simple basis set (delta function)
unbalance_events = false;
unbalance_classes = false;
whiten = false;
whiten_test_using_train = false;
resample_method = 'resample';
reproduce = false;
restrict_trainset = [];
for c=1:numel(methods)
    if any(strcmpi(methods{c},{'linear', 'quadratic', 'diagLinear', 'diagQuadratic', 'mahalanobis', 'svm'}))
        method = methods{c};
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
    if any(strcmpi(methods{c},{'induced','compute_induced'})) == 1
        compute_induced = true;
    end
    if any(strcmpi(methods{c},{'subtr_bin', 'subtr_indiv'})) == 1
        subtr_method = methods{c};
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
            if isnan(basis_sigma) || isempty(basis_sigma)
                basis_sigma = 1;
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
    if any(strcmpi(methods{c},{'resample', 'downsample', 'average_timebin'}))
        resample_method = methods{c};
    end
    if strncmpi(methods{c},'restrict_trainset',17)
        restrict_trainset = str2double(strrep(methods{c}(18:end),'-',','));
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
if compute_performance && (numel(unique([trainevents{:}])) < numel([trainevents{:}]) || numel(unique([testevents{:}])) < numel([testevents{:}]))
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
if ~compute_performance
    measuremethod = 'P(1stclass)'; % showing the average probability of the trial being classified as the first class in the class definition
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

% Determine bundle name and/or electrode selection
[channelset, bundlename_or_bundlelabels] = return_channel_bundle(channelset);

% a results folder for this electrode group
outpath = fullfile(outpath, channelset);

% a results folder for iterations / random permutations
if iterate && randomize_labels
    outpath = fullfile(outpath, 'randperm');
elseif iterate
    outpath = fullfile(outpath, 'iterations');
end

% create folder if it does not exist
if ~exist(outpath,'dir')
    mkdir(outpath);
end

% load data and determine output name
for cFile = 1:numel(filenames)
    msettings = [];
    msettings.channelpool = bundlename_or_bundlelabels;
    msettings.erp_baseline = erp_baseline{cFile};
    msettings.resample_eeg = resample_eeg; % NOTE: this line is different in TFR! Resampling here...
    msettings.resample_method = resample_method;
    msettings.do_csd = do_csd;
    msettings.clean_data = clean_muscle;
    msettings.clean_window = clean_window;
    msettings.shuffle_trials = true;
    [FT_EEG(cFile), filenames{cFile}, chanlocs{cFile}] = read_raw_data(filepath,filenames{cFile},outpath,msettings);
    % randomize labels for first level random permutation testing. NOTE: permuting all labels regardless of the conditions in the experiment
    if randomize_labels
        FT_EEG(cFile).trialinfo = FT_EEG(cFile).trialinfo(randperm(numel(FT_EEG(cFile).trialinfo)));
    end
end

% duplicate data for testing if only one file is available
if numel(filenames) == 1
    FT_EEG(2) = FT_EEG;
    chanlocs{2} = chanlocs{1};
end

% extract trialinfo, balance bins and restrict trial numbers
for cSet = 1:2
    thisCondSet = get_this_condset(condSet,cSet);
    if unbalance_events
        trialinfo{cSet} = FT_EEG(cSet).trialinfo; % use this for makefolds
    else
        % bin/balance dataset (default action, this is not to achieve actual binnning, it just applies within-class balancing of conditions)
        FT_EEG_BINNED(cSet) = compute_bins_on_FT_EEG(FT_EEG(cSet),thisCondSet,'trial','original');
        trialinfo{cSet} = FT_EEG_BINNED(cSet).trialinfo; % WE USE THIS BINNED TRIALINFO FOR MAKEFOLDS, THE ORIGINAL SETINDEX IS RECOVERED LATER ON!
        oldindex{cSet} = FT_EEG_BINNED(cSet).oldindex;   % OLDINDEX IS USED FOR RECOVERY TO ORIGINAL SETINDEX
        % a little hack to assign event number -99 to unbalanced events in FT_EEG (i.e. by selecting those events from thisCondSet that are not in FT_EEG_BINNED)
        origindex = find(ismember(FT_EEG(cSet).trialinfo,[thisCondSet{:}])); % all indices
        keepindex = [FT_EEG_BINNED(cSet).oldindex{:}]; % expand to get indices belonging to bins after balancing
        removeindex = setdiff(origindex,keepindex);
        FT_EEG(cSet).trialinfo(removeindex) = -99;
    end
    % limit trial numbers in each class of the training set if so desired, slightly complicated due to balancing
    if cSet == 1 && ~isempty(restrict_trainset)
        removeindex = [];
        nClasses = numel(thisCondSet);
        for cClass = 1:numel(thisCondSet)
            nCodesInClass = numel(thisCondSet{cClass});
            if unbalance_events % don't care about event codes
                restrictN = floor(restrict_trainset/nClasses);
            else % distribute evenly across instances in this class
                restrictN = floor((restrict_trainset/nClasses)/nCodesInClass); % trialindex is binned, so divide by number of event codes in each class
            end
            origindex = find(ismember(trialinfo{cSet},thisCondSet{cClass}));
            if numel(origindex) > restrictN
                removeindex = [removeindex origindex(restrictN+1:end)];
            else
                disp('WARNING: There are fewer trial instances in this class than the number you want to limit them by, so keeping the original number.');
            end
        end
        % assign event number -99 to remove events that exceed the restriction
        trialinfo{cSet}(removeindex) = -99;
        if unbalance_events
            FT_EEG(cSet).trialinfo = trialinfo{cSet}; % also put the new trialinfo back in FT_EEG
        else
            % a little hack to assign event number -99 to unbalanced events in FT_EEG (i.e. by selecting those events from thisCondSet that were not used to compute FT_EEG_BINNED)
            origindex = find(ismember(FT_EEG(cSet).trialinfo,[thisCondSet{:}])); % all indices before restricting
            keepindex = [ FT_EEG_BINNED.oldindex{ismember(FT_EEG_BINNED(cSet).trialinfo,[thisCondSet{:}])}]; % expand to get indices belonging to bins after restricting
            removeindex = setdiff(origindex,keepindex);
            FT_EEG(cSet).trialinfo(removeindex) = -99;
        end
    end
    % compute ERPs (baseline corrected, resampled, channels already selected)
    FT_ERP{cSet} = compute_erp_on_FT_EEG(FT_EEG(cSet),thisCondSet,'trial','bin');
    % keep track of channels and time line
    channels{cSet} = FT_EEG(cSet).label;
    times{cSet} = FT_EEG(cSet).time;
end

% if testing and training are on different files, check consistency
if numel(filenames) > 1 && ~all(strcmpi(channels{1},channels{2}))
    error('The electrodes do not occur in the same order in testing and training, some coding required to fix this...');
end

% Generate indices for folds to do training and testing
[setindex{1}, setindex{2}, nFolds] = make_folds(trialinfo{1},trialinfo{2},condSet,nFolds);

% remove any duplicate trial instances in the test set when only determining labels (given that
% the same class labels were probably entered for class 1 and class 2)
if ~compute_performance
    allindex = unique([setindex{2}{1}; setindex{2}{2}]);
    setindex{2}{1} = allindex(1:end/2);
    setindex{2}{2} = allindex(end/2:end);
end

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

% save ERPs
if ~iterate % do not save ERPs for iterations or every iteration would contain the ERPs
    fullfilename = [ outpath filesep filename ];
    save(fullfilename, 'FT_ERP','-v7.3');
end
clear FT_EEG_BINNED FT_ERP; % save memory by clearing

% between-class balancing: balance class instances by oversampling (only applied to training set)
if unbalance_classes
    wraptext('Between-class balancing is OFF. Make sure you know what you are doing, this can have undesirable effects on classifier bias when you have unevevenly represented stimulus classes in your design.',80);
else
    wraptext('Between-class balancing is ON, so that each stimulus class is evenly represented in the training set. If this is undesirable behavior, specify ''unbalance_classes'' in your methods.',80);
end

% within-class balancing: balance events within classes
if unbalance_events
    wraptext('Within-class balancing is OFF, such that an unequal distribution of event codes is allowed to contribute to each stimulus class. Make sure you know what you are doing, this can have undesirable effects on how you interpret your results.',80);
else % HERE WE RECOVER THE ORIGINAL SETINDEX!
    [setindex{1}, setindex{2}] = unpack_binned(setindex{1}, setindex{2}, oldindex{1}, oldindex{2}); % unpack setindex{1} and setindex{2} to get back the original index
    wraptext('Within-class balancing is ON, such that event codes are evenly represented within each stimulus class. If event codes are very unevenly represented in your data, this can result in the loss of many trials. It does however, enforce a balanced design, which is important for interpretation. If this is undesirable behavior, specify ''unbalance_events'' in your methods.',80);
end

% generate folds, compute some stuff on them and save them temporarily before running MVPA
settrialinfo = [];
settrialindex = [];
for cFld=1:nFolds
    % select trials from each dataset
    for cSet = 1:2
        % select trials belonging to this subset from original dataset
        trialindex = vertcat(setindex{cSet}{cFld,:});
        FT_EEG_2use = select_trials_from_FT_EEG(FT_EEG(cSet),trialindex);
        % get some goodies
        trialinfo{cSet} = FT_EEG_2use.trialinfo;
        origtrialindex{cSet} = FT_EEG(cSet).origindex(trialindex)';
        % get relevant condSet for computation
        condSet_2use = get_this_condset(condSet,cSet);
        % compute induced?
        if compute_induced
            if strcmpi(subtr_method,'subtr_bin')
                FT_EVOKED = compute_erp_on_FT_EEG(FT_EEG_2use,condSet_2use,'trial','bin');
            elseif strcmpi(subtr_method,'subtr_indiv')
                % subtracts the erp from each condition in a condSet
                FT_EVOKED = compute_erp_on_FT_EEG(FT_EEG_2use,condSet_2use,'trial','indiv');
            end
            if numel(trialindex) < (25*numel(condSet_2use)) % default action if there are not enough trials in the set to compute reliable ERPs
                disp('WARNING: fewer than 25 trials to compute an ERP for subtraction, using spline ERP for subtraction to be safe.');
                FT_EVOKED = compute_spline_on_FT_EEG(FT_EVOKED);
            end
            % subtract out the evoked responses
            FT_EEG_2use = subtract_evoked_from_FT_EEG(FT_EEG_2use,FT_EVOKED);
            clear FT_EVOKED;
        end
        % do binning if needed
        if (bintrain && cSet == 1) || (bintest && cSet == 2)
            FT_EEG_2use = compute_bins_on_FT_EEG(FT_EEG_2use,condSet_2use,'trial','original');
        end
        % get the goodies and get rid of overhead
        FT_EEG_2use = fix_dimord(FT_EEG_2use,'rpt_chan_time');
        % temporarily save fold data
        [~,tmpf,~] = fileparts(tempname); % generates random filename
        fnames{cFld,cSet} = [filepath filesep '..' filesep 'RAW_EEG_' filenames{1} '_fold' num2str(cFld) '_' num2str(cSet) '_' tmpf  '.mat'];
        save(fnames{cFld,cSet},'-v7.3','-struct','FT_EEG_2use');
        clear FT_EEG_2use;
    end
    % FYI, store for every fold
    settrialinfo = [settrialinfo; trialinfo];
    settrialindex = [settrialindex; origtrialindex];
end % end folds loop in which temp files are created

clear FT_EEG; % clear the dataset as we don't need it in memory during analyses

% and run models, with as little overhead as possible
for cFld=1:nFolds
    
    % FYI
    fprintf(1,['fold: ' num2str(cFld) '\n']);
    
    % load data, clear memory and delete obsolete files.
    train_FT_EEG = load(fnames{cFld,1}); % FT_EEG.trial = trial * channel * time
    test_FT_EEG = load(fnames{cFld,2});
    
    train_condSet = get_this_condset(condSet,1);
    test_condSet = get_this_condset(condSet,2);
    clear FT_EEG_foldsave;

    % fix dimord for whiten_FT_EEG and balance_FT_EEG functions, could be standardized to save memory peaks
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
    % NEED TO CHECK INSIDE balance_FT_EEG IF FT_EEG IS NOT ALREADY BINNED:
    if ~unbalance_classes && ~unbalance_events % only works when events within classes are balanced, otherwise just discard, not that important
        train_FT_EEG = balance_FT_EEG(train_FT_EEG,train_condSet,whiten);
    end
    
    % fix dimord for BDM_and_FEM_FT_EEG function, could be standardized to save memory peaks
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
    msettings.fname = fnames{cFld,1};
    
    % run BDM and FEM
    [BDM, FEM] = BDM_and_FEM_FT_EEG(train_FT_EEG,test_FT_EEG,train_condSet,test_condSet,msettings);
    
    % delete obsolete data and files
    clear train_FT_EEG test_FT_EEG train_condSet test_condSet
    delete(fnames{cFld,1});
    delete(fnames{cFld,2});

    % some BDM/FEM specific stuff
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
        FEM_ClassOverT(cFld,:,:) = tuning_from_matrix(FEM.C2_average,'slope',crossclass); % fld x t1 x t2 -> maybe move this to BDM_and_FEM_FT_EEG!
        FEM_WeightsOverT(cFld,:,:,:) = FEM.WeightsOverTime; % fld x time x elec x channel_response
        FEM_C2_average(cFld,:,:,:) = FEM.C2_average; % fld x time x time x channel_response
        FEM_C2_percondition(cFld,:,:,:,:) = FEM.C2_percondition; % fld x time x time x cond x channel_response
    end
end % end folds loop

% mean data, average over folds where applicable
BDM = [];
FEM = [];
if do_BDM
    % average over folds:
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
    % currently not possible to save confidence of FEM, maybe implement later
    clear FEM_*;
    BDM_CONF = [];
end

% save some settings so we know wtf just happened
settings.nconds = numel(condSet);
settings.nfolds = nFolds;
settings.filenames = filenames;
settings.crossclass = crossclass;
settings.erp_baseline = erp_baseline;
settings.resample_eeg = resample_eeg; % NOTE: this line is different in TFR! Resampling here...
settings.resample_method = resample_method;
settings.clean_window = clean_window;
settings.BDM = do_BDM;
settings.FEM = do_FEM;
settings.basis_set_sigma = basis_sigma;
settings.method_string = setmethod;
settings.channelset = channelset;
settings.channels = channels;
settings.chanlocs = chanlocs;
settings.times = times;
settings.frequency = 'none:_raw_eeg';
settings.dimord = 'time_time';
settings.measuremethod = measuremethod;
settings.compute_performance = compute_performance;
settings.save_confidence = save_confidence;
settings.restrict_trainset = restrict_trainset;
settings.trialinfo = settrialinfo;
settings.trialindex = settrialindex;
settings.condset = condSet;
settings.csd_transform = do_csd;
settings.bintrain = bintrain;
settings.bintest = bintest;
settings.unbalance_events = unbalance_events;
settings.unbalance_classes = unbalance_classes;
settings.whiten = whiten;
settings.whiten_test_using_train = whiten_test_using_train;

% SAVE RESULTS
% count filenames from 001 onwards if computing under permutation or iteration
if iterate
    fullfilename = find_filename(outpath,filename);
    save(fullfilename, 'FEM', 'BDM', 'BDM_CONF', 'settings', '-v7.3'); % not saving ERP to save memory
else
    fullfilename = fullfile(outpath, filename);
    save(fullfilename, 'FEM', 'BDM', 'BDM_CONF', 'settings', '-v7.3', '-append'); % this file also contains the ERPs, so append
end
warning('on','all');

function findfile = find_filename(path,filename)
c = 1;
findfile = fullfile(path, sprintf([filename '_PERM%04d'], c));
while numel(dir([findfile '.*']))>0
    c = c + 1;
    findfile = fullfile(path, sprintf([filename '_PERM%04d'], c));
end