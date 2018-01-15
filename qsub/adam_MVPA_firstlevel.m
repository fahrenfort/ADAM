function adam_MVPA_firstlevel(cfg)
% ADAM_MVPA_FIRSTLEVEL runs a first level MVPA analysis on single-subject raw EEG or MEG data.
% Output can be used as input for adam_compute_group_MVPA and adam_compute_group_ERP.
% 
% Use as:
%   adam_MVPA_firstlevel(cfg);
% 
% The function accepts as data formats: epoched raw/preprocessed files imported into EEGLAB (.set
% and .fdt) or Fieldtrip (.mat)
% 
% EEGLAB input files should be saved in EEGLAB native format (.set and .fdt) and contain an epoched
% EEGLAB data structure containing trials, channels and time. The events should contain triggers
% that uniquely refer to the conditions in the experiment. This format should be the default after
% running EEGLAB's pop_epoch(). ADAM extracts numerical trigger values that occur at timepoint 0 in
% the EEG.event.type field (so the triggers around which the trial was epoched). These are the
% trigger values that adam_MVPA_firstlevel uses for trial selection.
% 
% Fieldtrip input files should be saved in Matlab native format (.mat) and contain epoched Fieldtrip
% raw (see ft_datatype_raw) or timelock data (see ft_datatype_timelock), as long as they containing
% dimensions trial, time and channel, as well as a trialinfo array containing numerical trigger values that
% specify the condition of each trial. The fields that occur in a Fieldtrip data structure can
% either be saved as separate variables in the Matlab file (e.g. as happens when using the Matlab
% save command to save a data structure using the '-struct' keyword), or contain the struct variable
% itself (as when saving the data structure without the '-struct' keyword). In the latter case, the
% struct can have any name.
% 
% The cfg (configuration) variable serves to specify all required parameters for the MVPA
% analysis, including the directory where the input data is stored, the filenames (to run a batch of
% subjects in one go), and the directory where the output data should be saved.
% 
% The cfg should be specified as follows:
% 
%       cfg.datadir                = string specifiying the directory where the input files are
%                                    located; 
%       cfg.filenames              = N by 1 cell array containing N strings, in which each string
%                                    contains the base of a filename (e.g. cfg.filenames =
%                                    {'subject1' 'subject2'};). Do not add a file extension to the
%                                    file name, ADAM will automatically look for EEGLAB .set files,
%                                    if it cannot find these it will look for FieldTrip .mat files.
%                                    It is also possible to train your classifier on one input file
%                                    and test on a different input file (e.g. for cross-condition,
%                                    cross-task or cross-subject classification); this also requires
%                                    a N by 1 cell array, but now each string contains both
%                                    filenames (without extension) separated by a semi-colon (e.g.
%                                    cfg.filenames = {'subject1_1;subject1_2' 'subject2_1;subject2_2'};
%                                    in this example the classifier trains on the files that end in
%                                    _1  and tests on the files that end in _2)
%       cfg.class_spec             = N-class cell array, each cell containing a string with the
%                                    trigger values that are contained in that class, separated by
%                                    comma (e.g. cfg.class_spec{1} = '1,2,3'; cfg.class_spec{2} =
%                                    '4,5,6';). The function cond_string allows easy specification
%                                    of trigger values into the string format required for
%                                    adam_MVPA_firstlevel using variables containing integer arrays
%                                    (e.g. cfg.class_spec{1} = cond_string(faces); cfg.class_spec{2}
%                                    = cond_string(houses);). See cond_string for details. It is
%                                    also possible to define different triggers for training and
%                                    testing. In this case, train and test triggers need to be
%                                    separated by a semi-colon (e.g. cfg.class_spec{1} = '1,2;3,4';
%                                    cfg.class_spec{2} = '5,6;7,8'; in this example the classifier
%                                    trains on triggers (1,2) for class 1 and (5,6) for class 2,
%                                    while testing on (3,4) for class 1 and (7,8) for class 2).
%       cfg.balance_triggers       = 'yes' (default); balances triggers to achieve within-class
%                                    balancing, so that each class contains an equal amount of
%                                    trigger values (discarding leftover triggers); other option is
%                                    'no', which will use all triggers in the data, but be sure your
%                                    design is balanced! In general, we strongly recommend the default
%                                    'yes'. If desired, you can manually specify your own ratio of
%                                    each trial type in your class definition, by duplicating
%                                    trigger values: cfg.class_spec{1} = {'1,1,2,3'};
%                                    cfg.class_spec{2} = {'4,4,5,6'}.
%       cfg.balance_classes_method = 'oversample' (default; or 'undersample'); whether to
%                                    over/undersample classes in the training set to achieve
%                                    cross-class balancing; oversampling duplicates trials in the
%                                    training set.
%       cfg.class_type             = 'linear' (default); classifier type, e.g. 'linear' or
%                                    'diaglinear'; for other options see FITCDISCR (default Matlab
%                                    discriminant analysis function, which ADAM uses at its core)
%       cfg.class_method           = 'accuracy' (default); this the standard classification
%                                    metric which computes balanced accuracy (first computing
%                                    accuracy for each classs, and then averaging across classes);
%                                    other options are:  'hr-far' (hit rate minus false alarm
%                                    rate),'dprime' (d'),'hr' (hit rate),'far' (false alarm
%                                    rate),'mr' (miss rate),'cr' (correct rejection rate), in those
%                                    cases ADAM assumes that that the first class contains the
%                                    'signal' stimulus and the second class contains the 'noise'
%                                    stimulus. In the near future other measures will be
%                                    implemented, such as AUC and information prevalence.
%       cfg.crossclass             = 'no'; (default) whether ('yes') or not ('no') to apply
%                                    time-by-time cross-classification, yielding temporal
%                                    generalization matrices; specifying 'no' will simply compute
%                                    the diagonal of such a matrix (training and testing are done on
%                                    the same time points). Note that 'yes' will drastically
%                                    increase computation time, for exploratory analyses it may be
%                                    wise to down-sample prior to classification (see cfg.resample
%                                    below).
%       cfg.nfolds                 = integer specifying the number of folds for cross-validation
%                                    (default: 10);
%       cfg.model                  = 'BDM' (default); BDM performs a backward decoding model; FEM
%                                    performs a forward encoding model. One can also run both models
%                                    simultaneously by separating them with a comma (cfg.model =
%                                    'BDM,FEM';). It is only sensible to run a FEM when the classes
%                                    are thought to form a continuum, such as positions on a circle,
%                                    colors on a color wheel, or orientation.
%       cfg.raw_or_tfr             = 'raw' (default); you can either perform MVPA on raw EEG/MEG
%                                    data, or on time frequency representations of the data
%                                    (cfg.raw_or_tfr = 'tfr';) When 'tfr' is specified, ADAM first
%                                    peforms a time-frequency decomposition of the data (using
%                                    default settings), and subsquently peforms  decoding on the
%                                    single trial power values of each of these frequency bands. See
%                                    below for specification of frequency bands.
%       cfg.channels               = 'ALL_NOSELECTION'; (default) will use all channels/electrodes
%                                    in the data; note that if you still have EOG,ECG,etc. channels in
%                                    your data set, it will also include these, which may or may not
%                                    be desirable. If this is not desirable you can either remove
%                                    these channes from the input files prior to running ADAM, or
%                                    specify which channels/electrodes to include. There are two
%                                    ways of doing that. The first is using one of the default
%                                    electrode selections, e.g.: 'OCCIP', 'TEMPORAL', 'PARIET' or
%                                    'FRONTAL' (e.g. cfg.channels = 'OCCIP';). However, these will
%                                    only include electrodes from a standard 64-electrode BioSemi
%                                    10-20 system. It is easy to create your own selection by
%                                    modifying the select_channels function. The second method is to
%                                    directly specify the electrodes to include using a comma
%                                    separated list (e.g. cfg.channels = 'O1,O2,Iz,Oz';). type help
%                                    select_channels for more details.
%       cfg.clean_window           = [int int]; to remove muscle artefacts prior to decoding.
%                                    Specify a time window to inspect for muscle artefacts using
%                                    [begin,end]; always in SECONDS, e.g. cfg.clean_window =
%                                    [.25,1]; removes trials containing muscle artifacts between 250
%                                    and 1000 ms. For each subject, a .txt file and a .png graphic
%                                    will be saved in its results folder, showing a graphical
%                                    depiction of the artefacts that were removed (.png), as well as
%                                    a list of the trial numbers that were removed (.txt).
%       cfg.resample               = 'no' (default); or specify an integer to which to downsample
%                                    your data; this is recommended if you have a high sampling rate
%                                    (e.g. >250 Hz) and you want to do perform cross-classification
%                                    (temporal generalization, see under cfg.crossclass above). When
%                                    downsampling, it is recommended that the the original sampling
%                                    frequency is a multiple of the frequency to which to downsample
%                                    (e.g. if the original sampling frequency is 512Hz, do
%                                    cfg.resample = 256; or cfg.resample = 128; etc). When
%                                    time-frequency representations are computed, resampling is
%                                    applied only after time-frequency decomposition (so the TFRs
%                                    are computed on the original data).
%       cfg.erp_baseline           = 'no' (default); or specify a time window according to
%                                    [begin,end]; always in SECONDS, e.g. cfg.erp_baseline =
%                                    [-.25,0];
%       cfg.tfr_baseline           = 'no' (default); or specify a time window according to
%                                    [begin,end]; always in SECONDS, e.g. cfg.tfr_baseline =
%                                    [-.45,-.2];
%       cfg.frequencies            = '2:2:30' (default); takes frequencies from 2 to 30 Hz in steps
%                                    of 2; this should be a string, and only applies when
%                                    cfg.raw_or_tfr is set to 'tfr';
%       cfg.tfr_method             = 'total' (default); computes total power, alternatives are
%                                    'induced', which subtracts the erp from each trial (separately
%                                    performed on train and test data) and 'evoked'. When specifying
%                                    'evoked', ADAM averages groups of trials prior to computing
%                                    power. The averages are computed according to the class
%                                    definitions. E.g. if cfg.class_spec = '1,2,3'; it will compute
%                                    time frequency representations on the average of a triplet of
%                                    trials that have trigger codes 1, 2 and 3. One can increase the
%                                    number of trials that go into such an evoked average using the
%                                    class definition, e.g. when using cfg.class_spec{1} =
%                                    '1,1,1,2,2,2,3,3,3'; then each average on which power is
%                                    computed for class 1 will contain 9 rather than 3 trials.
%       cfg.bintrain               = 'no' (default); if 'yes', averages across triggers within a
%                                    class on the training side (for TFR data this only happens
%                                    after TF decomposition, do not use binning for 'evoked', which
%                                    averages prior to TF decomposition)
%       cfg.bintest                = 'no' (default); if 'yes', averages across triggers witin a
%                                    class on the testing side (for TFR data this only happens after
%                                    TF decomposition, do not use binning for 'evoked', which
%                                    averages prior to TF decomposition).
%       cfg.qsub                   = Struct that allows running these analyses on a Linux computing
%                                    cluster. When specifying this variable, the analysis will not
%                                    run locally but attempt to generate qsub files to be submitted
%                                    on a Linux computing cluster. The function will also generate
%                                    bash files to execute these qsub files using a single bash
%                                    command. In order to use this function, you need to compile the
%                                    classify_RAW_eeglab_data (for 'raw') and/or
%                                    classify_TFR_from_eeglab_data (for 'tfr') functions on the
%                                    linux cluster. These can be found in the /ADAM/eeg_mvpa folder.
%                                    The qsub struct can have a number of fields specifying the
%                                    parameters used when creating the qsub files, e.g.
%                                    qsub.walltime = '28:59:00'; 
%                                    qsub.cores = 16; 
%                                    qsub.maxcores = 10; 
%                                    qsub.mem = '64gb'; 
%                                    qsub.repeat = 1;
%                                    qsub.keep_together = true; 
%                                    qsub.qsubdir = '/Users/VU-MBP/lisa_remote'; 
%                                    qsub.functionpath = '$HOME/ADAM_master/eeg_mvpa';
%                                    qsub.qsubdir is the directory to which the linux server is
%                                    mapped on your local computer. It will assume the presence of a
%                                    '$HOME' text string in the cfg.datadir directory, which it
%                                    replaces with qsubdir when writing qsub files. This is useful
%                                    when you want to run qsub file generation locally, but write
%                                    the files to a remote server. qsub.functionpath specifies in
%                                    which directory on the linux cluster the compiled files can be
%                                    found (classify_RAW_eeglab_data for and/or
%                                    classify_TFR_from_eeglab_data).
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Example cfg with all default options used:
%
% cfg = [];
% cfg.datadir       = 'path/to/data/;
% cfg.filenames     = {'subject01'
%                      'subject02'
%                      'subject03'
%                      'subject04'
%                      'subject05'};
% cfg.class_spec{1} = '1,2,3'; % face-triggers
% cfg.class_spec{2} = '4,5,6'; % house-triggers
%
% adam_MVPA_firstlevel(cfg);
%
% --> This analysis performs backward decoding on raw EEG data using all scalp channels assuming
%     10-05 labels; classifier is trained on separating houses from faces; training and testing are
%     done on the same time points (i.e. the "diagonal", no temporal generalization /
%     cross-classification), using 10-fold cross-validation, with balanced triggers and oversampled
%     balanced classes. No baseline correction is done.
% 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Example cfg with some custom ingredients:
%
% cfg = [];
% cfg.datadir       = 'path/to/data/;
% cfg.filenames     = {'subject01_c1,subject01_c2'
%                      'subject02_c1,subject02_c2'
%                      'subject03_c1,subject03_c2'
%                      'subject04_c1,subject04_c2'
%                      'subject05_c1,subject05_c2'};
% houses = 1:3;
% faces = 4:6;
% masked = [1 4];
% unmasked = [2 5]
% scrambled = [3 6]
% cfg.class_spec{1} = cond_string(masked,faces);   % masked face-triggers
% cfg.class_spec{2} = cond_string(masked,houses);  % masked house-triggers
% cfg.cross_class   = 'yes';
% cfg.raw_or_tfr    = 'tfr';
% cfg.frequencies   = '8:2:12';
% cfg.tfr_method    = 'induced';
% cfg.tfr_baseline  = [.5 .2];
% cfg.channels      = 'OCCIPARIET';
%
% adam_MVPA_firstlevel(cfg);
%
% --> This analysis performs backward decoding on induced alpha power EEG data, only on 15
%     occipito-parietal channels; the classifier is trained on separating masked houses from masked
%     faces. Training on each time point is tested on every other time-point, yielding a temporal
%     generalization across time (GAT)  matrix. Note that this gives 3 GAT matrices, for each
%     frequency (8, 10 and 12). Induced power is baseline corrected using a pre-stimulus baseline
%     from -500 to -200 ms.
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Same cfg as above, only changing the following lines:
% 
% cfg.class_spec{1} = [ cond_string(masked,faces) ';' cond_string(unmasked,faces) ];
% cfg.class_spec{2} = [ cond_string(masked,houses) ';' cond_string(unmasked,houses) ];
% cfg.nfolds        = 1;
%
% adam_MVPA_firstlevel(cfg);
% 
% --> Same analysis as above, but now training on masked and testing on unmasked stimuli. Note that
%     cross-condition generalization does not require cross-validation (so setting nfolds to 1),
%     because a different data set is used for training and testing.
%
% part of the ADAM toolbox, by J.J.Fahrenfort, VU, 2017
% 
% See also POP_EPOCH, FT_DATATYPE_TIMELOCK, FT_DATATYPE_RAW, SELECT_CHANNELS, COND_STRING,
% FITCDISCR, ADAM_COMPUTE_GROUP_MVPA

% default values
channels = 'all';           % in 64-electrode BioSemi this uses all electrodes except the EOG electrodes, other options: 'ALL_NOSELECTION' for other aquisition systems or MEG, or for BioSemi 'OCCIP' only occipital, 'PARIET' only parietal etc, type help select_channels
nfolds = 10;                % trains on 90% (9/10) and tests on 10% (1/10)
crossclass = 'no';          % only trains and tests on the same time points
resample = 'no';            % does not resample the data
erp_baseline = 'no';        % [-.1,0] baselines from -100 to 0 ms
tfr_baseline = 'no';        % [-.5,-.1] baselines from -500 to -100 ms
frequencies = '2:2:30';     % '2:2:30' takes frequencies from 2 to 30 Hz in steps of 2
class_method = 'accuracy';  % computes classification accuracy (other options, e.g.: 'hr-far','dprime','hr','far','mr','cr', in those cases make sure that the first class is 'signal' and the second the 'noise')
class_type = 'linear';      % classifier type, e.g. 'linear' or 'diaglinear'
model = 'BDM';              % performs decoding rather than a forward encoding model
raw_or_tfr = 'raw';         % performs the analysis on the raw data rather than the time-frequeny data
balance_triggers = 'yes';   % balances triggers within each stimulus class, so that each class contains an equal amount of trigger values (discarding leftover triggers)
balance_classes_method = 'oversample'; % whether to oversample or undersample classes in the training set to achieve balancing
bintrain = 'no';            % average across triggers within a class on the training side
bintest = 'no';             % average across triggers witin a class on the testing side
savelabels = 'no';          % if 'yes', also saves the classifier labels
labelsonly = 'no';          % if 'yes', only saves the classifier labels (test set does not require labels in this case)
tfr_method = 'total';       % computes total power, alternative is 'induced' or 'evoked' ('induced' subtracts the erp from each trial, separately for train and test data, 'evoked' takes ERPs as input for TFR)
clean_window = [];

% unpack cfg
v2struct(cfg);

% do some checking
if ~exist('datadir','var')
    error('You need to specify in cfg.datadir where the data are located');
end
if ~exist('outputdir','var')
    error('You need to specify in cfg.outputdir where the results should be stored');
end
if ~exist('filenames','var')
    error('You need to specify a cell array in cfg.filenames containing the filenames containing the data of each subject');
end
if ~exist('class_spec','var')
    error('You need to specify in cfg.class_spec which trigger values go into which stimulus class used for training/testing');
end

% re-structure parameters to work with lower-level API
% settings string
if strcmpi(balance_triggers,'no')
    balance_triggers = 'unbalance_triggers';
else
    balance_triggers = '';
end
if strcmpi(balance_classes_method,'undersample')
    balance_classes_method = 'undersample';
elseif strcmpi(balance_classes_method,'oversample')
    balance_classes_method = 'oversample';
elseif strcmpi(balance_classes_method,'none')
    balance_classes_method = 'unbalance_classes';
end
if strcmpi(bintrain,'yes')
    bintrain = 'bintrain';
else
    bintrain = '';
end
if strcmpi(bintest,'yes')
    bintest = 'bintest';
else
    bintest = '';
end
if strcmpi(savelabels,'yes')
    savelabels = 'savelabels';
else
    savelabels = '';
end
if strcmpi(labelsonly,'yes')
    labelsonly = 'labelsonly';
else
    labelsonly = '';
end
if ~isempty(clean_window)
    clean_window = ['clean' num2str(clean_window)];
end
str_settings = cellarray2csvstring({class_method,class_type,model,balance_triggers,balance_classes_method,bintrain,bintest,tfr_method,savelabels,labelsonly,clean_window});
% other settings
if strcmpi(crossclass,'no') || isempty(crossclass)
    crossclass = '0';
else
    crossclass = '1';
end
if strcmpi(resample,'no') || isempty(resample)
    resample = '0';
elseif ~ischar(resample)
    resample = num2str(resample);
end
crossclass_resample = sprintf('%s,%s',crossclass,resample);
if strcmpi(erp_baseline,'no') || isempty(erp_baseline)
    erp_baseline = '0,0';
elseif ~ischar(erp_baseline)
    erp_baseline = sprintf('%f,%f',erp_baseline);
end
if strcmpi(tfr_baseline,'no') || isempty(tfr_baseline)
    tfr_baseline = '0,0';
elseif ~ischar(tfr_baseline)
    tfr_baseline = sprintf('%f,%f',tfr_baseline);
end
tfr_and_erp_baseline = sprintf('%s;%s',tfr_baseline,erp_baseline);
if isempty(frequencies)
    frequencies = '2:2:100';
end
if ischar(channels) && strcmpi(channels,'all')
    channels = 'ALL';
elseif iscell(channels)
    channels = cellarray2csvstring(channels);
end
if isempty(channels)
    channels = 'ALL_NOSELECTION';
end

% run analysis
if ~exist('qsub','var') %run local
    for cSubj = 1:numel(filenames)
        %try
        if strcmpi(raw_or_tfr,'raw')
            classify_RAW_eeglab_data(datadir,filenames{cSubj},outputdir,nfolds,channels,str_settings,crossclass_resample,erp_baseline,class_spec{:});
        elseif strcmpi(raw_or_tfr,'tfr')
            classify_TFR_from_eeglab_data(datadir,filenames{cSubj},outputdir,nfolds,channels,str_settings,crossclass_resample,tfr_and_erp_baseline,frequencies,class_spec{:});
        end
        %catch ME
        %    disp([ME.message ', skipping subject ' filenames{cSubj}]);
        %end
    end
else % or create qsub files
    if strcmpi(raw_or_tfr,'raw')
        create_qsub_files(qsub.functionpath,'classify_RAW_eeglab_data',qsub,datadir,filenames,outputdir,nfolds,channels,str_settings,crossclass_resample,erp_baseline,class_spec{:});
    elseif strcmpi(raw_or_tfr,'tfr')
        create_qsub_files(qsub.functionpath,'classify_TFR_from_eeglab_data',qsub,datadir,filenames,outputdir,nfolds,channels,str_settings,crossclass_resample,tfr_and_erp_baseline,frequencies,class_spec{:});
    end
end