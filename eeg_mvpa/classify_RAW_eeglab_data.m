function classify_RAW_eeglab_data(filepath,filenames,outpath,nFolds,channelset,method,crossclass_and_or_resample,erp_baselines,varargin)
% function classify_RAW_eeglab_data(filepath,filenames,outpath,nFolds,channelset,method,crossclass_and_or_resample,erp_baselines,varargin)
% Loads eeglab data (no time frequency decomposition) and performs a
% multivariate classification procedure.
% Filenames either contains a single filename for testing and training, or
% two filenames separated by a comma (the first for training, second for
% testing).
% outpath contains the folder where results should be stored (if empty
% defaults to filepath).
% nFolds is the number of cohorts in which the data is split up for
% the leave-one-out procedure (only when the same data are used for
% training and testing of course). When nFolds > 1, the algorithm divides
% trials over the folds such that each fold uses a fraction of
% (nFolds-1)/nFolds of all trials for training, and uses a fraction of
% 1/nFolds of all trials for testing. The training and testing procedure is
% applied to each fold separately, and the results of each of the folds is
% averaged at the end, ensuring that all data is tested once, without ever
% using the same data for training and testing. If some - but not all -
% trial labels are the same in the training and testing set
% (e.g. when some labels are used for training and testing, but other
% labels are unique to the training and/or testing set), the same
% procedure of splitting up trials for testing and training is followed
% separately for trials that are common to training and testing and the
% ones that are unique to training and testing, ensuring that their
% relative contribution during testing and training remains the same as it
% is in the original dataset. If you are using the same  dataset for
% training and testing, only set nFolds to 1 if training and testing labels
% are completely independent in the dataset. If there are separate sets for
% testing and training, nFolds defaults to 1.
% channelset determines the electrode subset that is used for testing. This
% can be done either numerically (assumes a 64 channel 10-20 system):
% 1 (uses all electrodes)
% 2 (uses occipital electrodes)
% 3 (uses parietal electrodes)
% 4 (uses frontal electrodes)
% 5 (uses temporal electrodes)
% 6 (uses occipitoparietal electrodes)
% 
% method specifies the method used for classification. Default is 'linear'
% for linear discriminant analysis. Other options are: 'diagLinear',
% 'mahalanobis' or 'quadratic'.
% The dependent classification measure can be specified in method,
% either using 'accuracy' (default), 'hr-far' or 'dprime' (last two only
% work when two categories are present, the first category is assumed to be
% the target present category, the second the target absent category. If
% more than two categories are present, the script defaults back to
% accuracy.
% Specify by adding to method, e.g. like this:
% method = 'linear,hr-far'
% If you want to output the individual labels that were assigned
% by the classifier to each of the trials rather than computing accuracy,
% specify 'labelsonly'. If no conditions are specified for the testing set,
% the algorithm goes through the entire testing set and outputs a label for
% each time point in each trial. No accuracy scores will be computed in
% this case. If conditions are specified for the test set, only conditions
% fitting this specification are labeled.
% In method you can also specify that the classifier should be trained
% under random permutation of the training labels by adding 'randperm',
% separated by a comma, like this: method = 'linear,randperm'. The
% results of computing under random permutation are stored in a folder
% called 'randperm' and counted so that they can be used for non-parametric
% hypothesis testing.
% Similarly, you can also run another complete iteration of the same data
% (all folds) without randomly permuting the labels, by specificying
% method = 'iterate'. The results of computing another iteration are stored
% in a folder called 'iterations' and counted so that they can subsequently
% be used for averaging to get a cleaner result.
% You can specify the number of times you want to execute a permutation or
% iteration through the create_qsub_files function, by setting the number
% of repeats using settings.repeat as an argument to that function.
% In method, you can also specify whether the trial labels that are
% specified for each condition should be binned (averaged) to generate new
% trials, or not. This is done by adding 'bin' to the method,
% separated by comma's, like this:
% method = 'linear,bin';
% If left out, no binning is applied (default). One can also specify
% 'bintrain' or 'bintest' if binning should only take place on the training
% or on the testing side.
% Other options for method are to specify whether to compute the Scalp
% Current Density (SCD) prior to classification. This is done by either
% specifying 'csd' or 'scd' (e.g. method = 'linear,csd') using the finite
% method from ft_scalpcurrentdensity from the fieldtrip toolbox.
% crossclass_and_or_resample specifies whether a cross-classification is
% performed over time, in which each time point is used to classify all
% other time points. This is very time consuming, better set to 0 for first
% analysis, but set to 1 if you want to get a better feel of the data. You
% can speed up cross classification by downsampling the data. This is done
% by specifying the new sampling rate, separating the value for whether or
% not to cross classify from the new sampling rate using a semicolon, like
% this:
% crossclass_and_or_resample = '1;128'; meaning that cross classification
% will be performed and that all data will be resampled to 128 Hz. If no
% value is specified to resample to, the data is left at the original
% sampling rate.
% erp_baselines allows you to perform baseline correction on the data.
% Specify in seconds [begin,end]. If specified as 0 or [0,0], no baseline
% correction will be applied (default). It is also possible to specify
% separate baselines for two input sets when separating them using a
% semicolon, like this: erp_baselines = '-.25,0; 0,0', which would apply
% a baseline correction of -250 to 0 ms only to the training set only.
% varargin is a variable set containing the trigger codes to select
% conditions from, specified as string with comma separated values:
% cond1 = '1,2,3';
% cond2 = '4,5,6';
% In that setup, all trials containing either a trigger code 1, 2 or 3 will
% be taken as belonging to condition 1, while all trials containing a 4, 5
% and 6 will be taken as belonging to condition 2.
% If different trigger codes are used in training versus testing, separate
% the condition specification by a semicolon:
% cond1 = '1,2;8,9';
% cond2 = '3;10';
% here 1, 2 (the first) would be for training and 8,9 (the second for
% testing) condition 1, while 3 would be training condition 2 and 10 for
% testing condition 2.
% Usage examples:
% classify_RAW_data('/home/eeg/TFR','subject1data1,subject1data2','home/eeg/mvpa/TFR',1,1,'linear,bintrials,hr-far',0,'-.2,0;-.4,-.2','1,2,3','4,5,6')
% classify_RAW_data('/home/eeg/TFR','subject1data','home/eeg/mvpa/TFR',4,2,'linear',0,0,'1;3','2;4')
% J.J.Fahrenfort, VU 2014,2015

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
measuremethod = 'accuracy';
randomize_labels = false;
iterate = false;
bintest = false;
bintrain = false;
labelsonly = false;
do_csd = false;
use_splines = false;
compute_induced = false;
do_FEM = false;
do_BDM = false;
clean_muscle = false;
clean_window = [];
save_labels = false;
basis_sigma = 1; % default width of basis set if not a delta, if this is empty, do a simple basis set (delta function)
unbalance_triggers = false;
unbalance_classes = false;
balance_classes_method = 'oversample';
detrend_eeg = false;
for c=1:numel(methods)
    if any(strcmpi(methods{c},{'linear', 'quadratic', 'diagLinear', 'diagQuadratic', 'mahalanobis'}))
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
        randomize_labels = true;
        % create a folder for random permutation
        outpath = [outpath filesep 'randperm'];
        mkdir(outpath);
    end
    if strcmpi(methods{c},'iterate')
        iterate = true;
        % create a folder for iterations
        outpath = [outpath filesep 'iterations'];
        mkdir(outpath);
    end
    if any(strcmpi(methods{c},{'induced','compute_induced'})) == 1
        compute_induced = true;
    end
    if any(strcmpi(methods{c},{'subtr_bin', 'subtr_indiv'})) == 1
        subtr_method = methods{c};
    end
    if any(strcmpi(methods{c},{'hr-far','dprime','hr','far','mr','cr'}))
        measuremethod = methods{c};
        if numel(condSet) ~= 2
            disp('Number of stimulus classes is unequal to 2, defaulting back to computing accuracy rather than hr-far');
            measuremethod = 'accuracy'; % defaulting back to accuracy
        else
            disp('computing sdt measure, assuming the first condition is target (signal) and second is non-target (noise)');
        end
    end
    if any(strcmpi(methods{c},{'labelsonly','onlylabels','labels_only'}))
        labelsonly = true;
        save_labels = true;
        crossclass = false; % cross class is useless when keeping all labels, we will run out of memory (time * time * trials)
        do_FEM = false;
        do_BDM = true;
        disp('Labels only works with BDM only');
    end
    if any(strcmpi(methods{c},{'csd','scd'}))
        do_csd = true;
    end
    if any(strcmpi(methods{c},{'splines','spline'}))
        use_splines = true;
    end
    if any(strcmpi(methods{c},{'FEM','do_FEM','forward'})) && ~labelsonly
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
    if any(strcmpi(methods{c},{'save_labels','savelabels'}))
        save_labels = true;
    end
    if any(strcmpi(methods{c},{'unbalance_triggers','unbalance','unbalanced'}))
        unbalance_triggers = true;
    end
    if any(strcmpi(methods{c},{'unbalance_classes'}))
        unbalance_classes = true;
        balance_classes_method = 'none';
    end
    if any(strcmpi(methods{c},{'undersample'}))
        unbalance_classes = false;
        balance_classes_method = 'undersample';
    end
    if any(strcmpi(methods{c},{'oversample'}))
        unbalance_classes = false;
        balance_classes_method = 'oversample';
    end
    if any(strcmpi(methods{c},{'detrend','detrend_eeg'}))
        detrend_eeg = true;
    end
end
if ~do_FEM && ~do_BDM
    do_BDM = true;
end

% check nFolds and condSet
if numel(filenames) > 1 && nFolds > 1
    disp('WARNING: You specified different filenames for training and testing, with more than 1 fold. Leave-one-out is not applicable here. Defaulting nFolds to 1.');
    nFolds = 1; % if you do want to cut up independent sets into folds, you can do so by turning this safety check off
end
if isempty(condSet)
    error('Cannot find usable trigger specification.');
end
% double condset for test and train if only one is specified
if size(condSet{1},1) == 1
    condSet = put_this_condset(condSet,condSet,2);
end
% if using same triggers for training and testing, increase nFolds
if numel(filenames) == 1 && nFolds == 1
    for cCondSet = 1:numel(condSet)
        if any(ismember(condSet{cCondSet}(1,:),condSet{cCondSet}(2,:))) && nFolds == 1
            nFolds = 10;
            wraptext('WARNING: You dirty double dipper! You are using the same data for testing and training without a leave-one-out procedure. Defaulting nFolds to 10 for crossvalidation.',80);
        end
    end
end
% check if condsets are not the same but overlapping, if so unbalance_triggers
for cCondSet = 1:numel(condSet)
    if ~all(ismember(condSet{cCondSet}(1,:),condSet{cCondSet}(2,:))) && any(ismember(condSet{cCondSet}(1,:),condSet{cCondSet}(2,:))) && numel(filenames) == 1 && ~unbalance_triggers
        unbalance_triggers = true;
        wraptext('WARNING: Some stimulus triggers overlap between train and test, overriding balance triggers option');
    end
end

% display stimulus classes
wraptext('These are the stimulus classes. Each row contains the trigger codes that go into a single class (first row training, second row testing):',80);
celldisp(condSet,'stimclass');

% Determine bundle name and/or electrode selection
[channelset, bundlename_or_bundlelabels] = return_channel_bundle(channelset);

% create a folder for this electrode group
outpath = [outpath filesep channelset];
if ~exist(outpath,'dir')
    mkdir(outpath);
end

% load data and determine output name
for cFile = 1:numel(filenames)
    % NOTE: this line is different in TFR! Resampling here...
    [FT_EEG(cFile), filenames{cFile}, chanlocs{cFile}]= read_raw_data(filepath,filenames{cFile},outpath,bundlename_or_bundlelabels,erp_baseline{cFile},resample_eeg,do_csd,clean_muscle,clean_window,true,detrend_eeg);
end
 % duplicate data for testing if only one file is available
if numel(filenames) == 1
    FT_EEG(2) = FT_EEG;
end

% extract some relevant trial info, training and testing data
for cSet = 1:2
    thisCondSet = get_this_condset(condSet,cSet);
    if unbalance_triggers
        trialinfo{cSet} = FT_EEG(cSet).trialinfo;
    else
        % bin/balance dataset (default action, this is not to achieve actual binnning, it just balances the dataset in case separate conditions still exist in each stimulus class)
        FT_EEG_BINNED(cSet) = compute_bins_on_FT_EEG(FT_EEG(cSet),thisCondSet,'trial','original');
        trialinfo{cSet} = FT_EEG_BINNED(cSet).trialinfo;
        oldindex{cSet} = FT_EEG_BINNED(cSet).oldindex;
        % remove non-balanced trials for ERP calculation
        tempbool = ones(size(FT_EEG(cSet).trialinfo)); tempbool([oldindex{cSet}{:}]) = 0; 
        FT_EEG(cSet).trialinfo(logical(tempbool)) = NaN;
        clear tempbool;
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

% randomize labels if desired. NOTE: permuting all labels of the input set,
% regardless of whether there is a separate testing set or not
if randomize_labels
    trialinfo{1} = trialinfo{1}(randperm(numel(trialinfo{1})));
end

% Generate indices for folds to do training and testing
[setindex{1}, setindex{2}, nFolds] = make_folds(trialinfo{1},trialinfo{2},condSet,nFolds,labelsonly);

% create file name based on actual folds and save ERPs
if numel(filenames) == 1
    filename = ['CLASS_PERF_' filenames{1} '_' num2str(nFolds) 'fold'];
else
    filename = ['CLASS_PERF_' filenames{1} '_' filenames{2}];
end
if ~(randomize_labels || iterate)
    fullfilename = [ outpath filesep filename ];
    save(fullfilename, 'FT_ERP','-v7.3');
end
clear FT_EEG_BINNED FT_ERP; % save memory by clearing

% balance class instances by oversampling or undersampling (only applied to training set)
if unbalance_classes
    wraptext('Please realize that stimulus classes are now UNBALANCED. Make sure you know what you are doing, this can have undesirable effects on classifier bias when you have unevevenly represented stimulus classes in your design.',80);
else
    % duplicate or eliminate stimulus classes from the training set
    for cFld=1:nFolds
        nEachClass = cellfun(@numel, setindex{1}(cFld,:));
        if strcmpi(balance_classes_method,'oversample')
            maxN = max(nEachClass);
            disp('Balancing classes by oversampling: duplicating class instances in the training set (default).');
        elseif strcmpi(balance_classes_method,'undersample')
            maxN = min(nEachClass);
            disp('Balancing classes by undersampling: eleminating class instances from the training set.');
        end
        for cClass = 1:numel(nEachClass);
            elements = setindex{1}{cFld,cClass};
            elements = repmat(elements,ceil(maxN/numel(elements)),1);
            setindex{1}{cFld,cClass} = elements(1:maxN);
        end
    end
    wraptext('Stimulus classes are now BALANCED by design, so that each stimulus class is evenly represented in the training set. If this is undesirable behavior, specify ''unbalance_classes'' in your methods.',80);
end

% balance triggers within classes
if unbalance_triggers
    wraptext('Please realize that triggercodes in a class are now UNBALANCED, such that an unequal distribution of triggercodes is allowed to contribute to each stimulus class. Make sure you know what you are doing, this can have undesirable effects on how you interpret your results.',80);
else
    % unpack setindex{1} and setindex{2} to get back the original index
    [setindex{1}, setindex{2}] = unpack_binned(setindex{1}, setindex{2}, oldindex{1}, oldindex{2});
    wraptext('Triggercodes in a class are now BALANCED by design, such that triggercodes are evenly represented within each stimulus class. If triggercodes are very unevenly represented in your data, this can result in the loss of many trials. It does however, enforce a balanced design, which is important for interpretation. If this is undesirable behavior, specify ''unbalance_triggers'' in your methods.',80);
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
            if numel(trialindex) < (25*numel(condSet_2use)) % default action if there are not enough trials in test
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
        % compute splines if specified
        if use_splines
            FT_EEG_2use = compute_spline_on_FT_EEG(FT_EEG_2use);
        end
        % get the goodies and get rid of overhead
        FT_EEG_2use = fix_dimord(FT_EEG_2use);
        alldata{cSet} = FT_EEG_2use.trial;
        labels{cSet} = make_group_labels(FT_EEG_2use.trialinfo, get_this_condset(condSet,cSet));
        clear FT_EEG_2use;
    end
    % FYI, store for every fold
    settrialinfo = [settrialinfo; trialinfo];
    settrialindex = [settrialindex; origtrialindex];
    
    % temporarily save folded data
    [~,tmpf,~] = fileparts(tempname); % generates random filename
    tempfilename{cFld} = [filepath filesep '..' filesep 'RAW_EEG_' filenames{1} '_fold' num2str(cFld) '_' tmpf  '.mat'];
    save(tempfilename{cFld},'alldata','labels','-v7.3');
    
end % end folds loop in which temp files are created

clear FT_EEG; % clear the dataset as we don't need it in memory during analyses

% and run models, with as little overhead as possible
for cFld=1:nFolds
    fprintf(1,['fold: ' num2str(cFld) '\n']);
 
    % settings for backward and/or forward modelling
    msettings.crossclass = crossclass;
    msettings.method = method;
    msettings.labelsonly = labelsonly;
    msettings.doBDM = do_BDM;
    msettings.doFEM = do_FEM;
    msettings.basis_sigma = basis_sigma;
    % load data, run analysis, clear memory and delete obsolete files
    % note that Matlab does not copy a matrix that is passed into a  function when that matrix 
    % is not modified inside that function. Rather, it creates a pointer to save memory :-)
    load(tempfilename{cFld}); % alldata = elec * time * trial
    [BDM, FEM] = EEG_backward_and_forward_matrix(alldata{2},alldata{1},labels{2},labels{1},msettings);
    % delete obsolete data and files
    clear alldata labels;
    delete(tempfilename{cFld});

    % and BDM/FEM specific stuff
    if do_BDM
        if save_labels && ~labelsonly
            BDM_labelMatrixOverT(cFld,:,:,:,:) = BDM.LabelsOverTime; % fld x t1 x t2 x response_matrix
        elseif save_labels && labelsonly
            BDM_labelMatrixOverT{cFld} = BDM.LabelsOverTime; % fld x trial x t1 (assigned_labels when method = 'labelsonly', nr of trials can be different for each fold)
        end
        if labelsonly
            BDM_ClassOverT(cFld) = NaN;
        else
            BDM_ClassOverT(cFld,:,:) = class_accuracy_from_matrix(BDM.LabelsOverTime,measuremethod,crossclass); % fld x t1 x t2
        end
        BDM_WeightsOverT(cFld,:,:) = BDM.WeightsOverTime; % fld x time x elec
        BDM_covPatternsOverT(cFld,:,:) = BDM.covPatternsOverTime; % fld x time x elec
        BDM_corPatternsOverT(cFld,:,:) = BDM.corPatternsOverTime; % fld x time x elec
    end
    if do_FEM
        if save_labels
            FEM_labelMatrixOverT(cFld,:,:,:,:) = FEM.LabelsOverTime;
        end
        % FEM_ClassOverT(cFld,:,:) = class_accuracy_from_matrix(FEM.LabelsOverTime,measuremethod); % fld x t1 x t2
        FEM_ClassOverT(cFld,:,:) = tuning_from_matrix(FEM.C2_average,'slope',crossclass); % fld x t1 x t2
        FEM_WeightsOverT(cFld,:,:,:) = FEM.WeightsOverTime; % fld x time x elec x channel_response
        FEM_C2_average(cFld,:,:,:) = FEM.C2_average; % fld x time x time x channel_response
        FEM_C2_percondition(cFld,:,:,:,:) = FEM.C2_percondition; % fld x time x time x cond x channel_response
    end
    
end % end folds loop

% mean data, average over folds where applicable
BDM = [];
FEM = [];
BDMLabelsOverTime = [];
FEMLabelsOverTime = [];
if do_BDM
    BDM.ClassOverTime = squeeze(mean(BDM_ClassOverT,1)); % t1 x t2 or trial x t1
    BDM.WeightsOverTime = squeeze(mean(BDM_WeightsOverT,1)); % time x elec x channel_response
    BDM.covPatternsOverTime = squeeze(mean(BDM_covPatternsOverT,1)); % time x elec
    BDM.corPatternsOverTime = squeeze(mean(BDM_corPatternsOverT,1)); % time x elec
    if save_labels
        BDMLabelsOverTime = squeeze(BDM_labelMatrixOverT); % fld x t1 x t2 x response_matrix OR fold x trial x t1 (assigned_labels when method = 'labelsonly')
    end
    clear BDM_*;
end
if do_FEM
    FEM.ClassOverTime = squeeze(mean(FEM_ClassOverT,1)); 
    FEM.WeightsOverTime = squeeze(mean(FEM_WeightsOverT,1));
    FEM.C2_average = squeeze(mean(FEM_C2_average,1));
    FEM.C2_percondition = squeeze(mean(FEM_C2_percondition,1));
    if save_labels
        FEMLabelsOverTime = squeeze(FEM_labelMatrixOverT);
    end
    clear FEM_*;
end

% save some settings so we know wtf just happened
settings.nconds = numel(condSet);
settings.nfolds = nFolds;
settings.filenames = filenames;
settings.crossclass = crossclass;
settings.erp_baseline = erp_baseline;
settings.clean_window = clean_window;
settings.detrend_eeg = detrend_eeg;
settings.BDM = do_BDM;
settings.FEM = do_FEM;
settings.basis_set_sigma = basis_sigma;
settings.method_string = setmethod;
settings.use_splines = use_splines;
settings.channelset = channelset;
settings.channels = channels;
settings.chanlocs = chanlocs;
settings.times = times;
settings.frequency = 'none:_raw_eeg';
settings.dimord = 'time_time';
settings.measuremethod = measuremethod;
settings.trialinfo = settrialinfo;
settings.trialindex = settrialindex;
settings.condset = condSet;
settings.csd_transform = do_csd;
settings.bintrain = bintrain;
settings.bintest = bintest;
settings.unbalance_triggers = unbalance_triggers;
settings.unbalance_classes = unbalance_classes;
settings.balance_classes_method = balance_classes_method;

% count filenames from 001 onwards if computing under permutation or iteration
if randomize_labels || iterate
    fullfilename = find_filename(outpath,filename);
    save(fullfilename, 'FEM', 'BDM', 'settings', '-v7.3');
else
    fullfilename = [ outpath filesep filename ];
    save(fullfilename, 'FEM', 'BDM', 'settings', '-v7.3', '-append'); % this file also contains the ERPs, so append
end
if save_labels
    if labelsonly
        save_var_under_different_name(fullfilename,BDMLabelsOverTime, 'BDM_LabelsOverTime', FEMLabelsOverTime, 'FEM_LabelsOverTime');
    else
        save_var_under_different_name(fullfilename,BDMLabelsOverTime, 'BDM_ConfusionMatrixOverTime', FEMLabelsOverTime, 'FEM_ConfusionMatrixOverTime');
    end
end
warning('on','all')

function fullfile = find_filename(path,filename)
c = 1;
fullfile = sprintf([path filesep filename '_PERM%03d'], c);
while numel(dir([fullfile '.*']))>0
    c = c + 1;
    fullfile = sprintf([path filesep filename '_PERM%03d'], c);
end