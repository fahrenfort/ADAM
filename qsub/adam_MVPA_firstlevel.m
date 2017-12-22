function adam_MVPA_firstlevel(cfg)
% ADAM_MVPA_FIRSTLEVEL runs a first level MVPA analysis on single-subject raw EEG or MEG data.
% Output can be used as input for ADAM_COMPUTE_GROUP_MVPA.
% 
% Use as:
%   adam_MVPA_firstlevel(cfg);
% 
% The function accepts as data formats: raw/preprocessed files imported into EEGLAB (.set and .fdt)
% or Fieldtrip (.mat)
% 
% EEGLAB files should have a EEG.data structure with dimensions: channels by time by trials, and an
% EEG.event structure containing triggers that uniquely refer to the different classes on which the
% MVPA analysis is performed (see below). This format should be the default after running EEGLAB's
% pop_epoch().
% 
% Fieldtrip files should have a data.trial cell array containing all trials, where each trial-cell
% has dimensions: channel by time, and a data.trialinfo array containing triggers that uniquely
% refer to the different classes on which the MVPA analysis is performed (see below). This format
% should be the default after running Fieldtrip?s ft_preprocessing() with a proper trial definition
% (see ft_trialfun).
% 
% The cfg (configuration) structure serves as input and contains all ingredients for the MVPA
% analysis, as well as the directory where the input data is stored, the filenames (to run a batch
% of subjects in one go), and the directory where the output data should be saved.
% 
% The cfg should be specified as follows:
% 
%       cfg.datadir                = string specifiying a directory where the raw data are located;
%       cfg.filenames              = a N by 1 cell array with N filenames, without .mat or .set
%                                    extension; It is also possible to train your classifier on one
%                                    input file and test on another input file (cross-condition or
%                                    cross-subject classification); this still requires a N by 1
%                                    cell array, where each entry is a string with both filenames
%                                    (without extension) separated by a semi-colon (see below for an
%                                    example)
%       cfg.class_spec             = a N-class cell array, each cell containing a string with the
%                                    trigger values corresponding to that class, separated by comma
%                                    (e.g. cfg.class_spec{1} = '1,2,3'; cfg.class_spec{2} =
%                                    '4,5,6';). Pay attention to the balancing of your stimulus
%                                    classes. The function COND_STRING combines integer arrays of
%                                    triggers, according to possible multiple levels of your
%                                    experimental design (e.g. faces/houses and
%                                    conscious/unconscious), into the cell-array format needed for
%                                    ADAM_MVPA_FIRSTLEVEL.
%       cfg.balance_triggers       = 'yes' (default); balances triggers to achieve within-class
%                                    balancing, so that each class contains an equal amount of
%                                    trigger values (discarding leftover triggers); other option is
%                                    'no', which will use all triggers in the data, but be sure your
%                                    design is balanced! In general, we highly recommend the default
%                                    'yes'. If desired, you can manually specify your own ratio of
%                                    each trial type in your class definition, by duplicating
%                                    trigger values: {'1,1,2,3'} {'4,4,5,6'}.
%       cfg.balance_classes_method = 'oversample' (default; or 'undersample'); whether to
%                                    over/undersample classes in the training set to achieve
%                                    cross-class balancing; oversampling results in some duplicate
%                                    trials.
%       cfg.class_type             = 'linear' (default); classifier type, e.g. 'linear' or
%                                    'diaglinear'; for other options see FITCDISCR (default Matlab
%                                    discriminant analysis function, which ADAM uses at its core)
%       cfg.class_method           = 'accuracy' (default); this the "standard" classification
%                                    metric; other options are:
%                                    ?hr-far','dprime','hr','far','mr','cr', in those cases make
%                                    sure that the first class is 'signal' and the second the
%                                    'noise'
%       cfg.crossclass             = 'no'; (default) whether ('yes') or not ('no') to apply
%                                    time-by-time cross-classification, yielding temporal
%                                    generalization matrices; specifying 'no' will simply compute
%                                    the diagonal of such a matrix (training and testing are done on
%                                    the same time points). Note that 'yes' will drastically
%                                    increase computation time!
%       cfg.nfolds                 = integer specifying the number of folds for cross-validation
%                                    (default: 10);
%       cfg.model                  = 'BDM' (default); BDM performs a backward decoding model; FEM
%                                    performs a forward encoding model.
%       cfg.raw_or_tfr             = 'raw' (default); you can either perform MVPA on raw M/EEG data,
%                                    or on single trial power of multiple frequency bands (option
%                                    'tfr'; see below for specification of frequency bands).
%       cfg.channels               = 'ALL_NOSELECTION'; specifying this will simply use all channels
%                                    available; note that if you still have EOG,ECG,etc. channels in
%                                    your data set, it will also include these, which may or may not
%                                    be desirable. Other options exist, which currently only hold
%                                    for data with a channel structure that contains labels
%                                    according to the 10-05 labeling convention: 'all' (default; all
%                                    scalp-channels, so excluding EOG,etc.), 'OCCIP' (only
%                                    occipital), 'PARIET' (only parietal) etc; type help
%                                    select_channels.
%       cfg.resample               = 'no' (default); or specify an integer to downsample your data;
%                                    this is especially recommended if you have a high sampling rate
%                                    (e.g. >500 Hz) and you want to do perform cross-classification
%       cfg.erp_baseline           = 'no' (default); or specify a time window according to
%                                    [begin,end]; always in SECONDS.
%       cfg.tfr_baseline           = 'no' (default); or specify a time window according to
%                                    [begin,end]; always in SECONDS.
%       cfg.frequencies            = '2:2:30' (default); takes frequencies from 2 to 30 Hz in steps
%                                    of 2; this should be a string, and only applies when
%                                    cfg.raw_or_tfr is set to 'tfr';
%       cfg.tfr_method             = 'total' (default); computes total power, alternative is
%                                    'induced' or 'evoked' ('induced' subtracts the erp from each
%                                    trial, separately for train and test data, 'evoked' takes ERPs
%                                    as input for TFR).
%       cfg.bintrain               = 'no' (default); if 'yes', averages across triggers within a
%                                    class on the training side
%       cfg.bintest                = 'no' (default); if 'yes', averages across triggers witin a
%                                    class on the testing side
%       cfg.savelabels             = 'no' (default); if 'yes', also saves the classifier labels
%       cfg.labelsonly             = 'no' (default); if 'yes', only saves the classifier labels
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Example cfg with all default options used:
%
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
%     done on the same time points (i.e. the "diagonal", no temporal generalization
%     cross-classification), using 10-fold cross-validation, with balanced triggers and oversampled
%     balanced classes. No baseline correction is done.
% 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Example cfg with some custom ingredients:
%
% cfg.datadir       = 'path/to/data/;
% cfg.filenames     = {'subject01_c1,subject01_c2'
%                      'subject02_c1,subject02_c2'
%                      'subject03_c1,subject03_c2'
%                      'subject04_c1,subject04_c2'
%                      'subject05_c1,subject05_c2'};
% cfg.class_spec{1} = '1,2,3'; % face-triggers
% cfg.class_spec{2} = '4,5,6'; % house-triggers
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
%     parieto-occipital channels; the classifier is trained on separating houses from faces in one
%     experimental condition (say, not-masked), and tested on another experimental condition (say,
%     masked). Note that cross-condition generalization does not require cross-validation, because a
%     different data set is used for training and testing. Training on each time point is tested on
%     every other time-point, yielding a temporal generalization across time (GAT)  matrix. Note 
%     that this gives 3 GAT matrices, for each frequency (8, 10 and 12). Induced power is baseline
%     corrected using a pre-stimulus baseline from -500 to -200 ms.
%
% part of the ADAM toolbox, by J.J.Fahrenfort, VU, 2017
% 
% See also ADAM_COMPUTE_GROUP_MVPA, POP_EPOCH, FT_PREPROCESSING, SELECT_CHANNELS, COND_STRING,
% FITCDISCR

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

% unpack cfg
v2struct(cfg);

% do some checking
if ~exist('datadir','var')
    error('You need to specify the directory where the data are located');
end
if ~exist('outputdir','var')
    error('You need to specify the directory where the results should be stored');
end
if ~exist('filenames','var')
    error('You need to specify a cell array containing the filenames containing the data of each subject');
end
if ~exist('class_spec','var')
    error('You need to specify the trigger values that go into each stimulus class used for training/testing');
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
str_settings = cellarray2csvstring({class_method,class_type,model,balance_triggers,balance_classes_method,bintrain,bintest,tfr_method,savelabels,labelsonly});
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