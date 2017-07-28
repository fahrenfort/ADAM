function adam_MVPA_firstlevel(cfg)
% This function runs a first level analysis based on the parameters
% specificied in the cfg struct. Unspecified parameters take on default
% values. Amongst others, the following core parameters can be specified:
% cfg.datadir = '/path/to/the/data';
% cfg.outputdir = '/path/to/the/results';
% cfg.filenames = {'subject1', 'subject2);           do not use extensions, just base of filename!
%                                                    for separate train and test files do: 
% cfg.filenames = {'subject1train;subject1test', 'subject2train;subject2test'};
% cfg.class_spec = {'1,2,3','4,5,6'};                triggers 1,2,3 are class 1, and 4,5,6 are class 2
%                                                    for separate train and test triggers, do:
% cfg.class_spec = {'1,2,3;11,12','4,5,6;13,14'};    trains on triggers 1,2,3 (class1) and 4,5,6 (class2),
%                                                    but tests on triggers '11,12' (class1) and '13,14' (class2)
%
% part of the ADAM toolbox, by J.J.Fahrenfort, VU, 2017

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
tfr_method = 'total';       % computes total power, alternative is 'induced' (subtracts the erp from each trial, separately for train and test data)

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
str_settings = sprintf('%s,%s,%s,%s,%s,%s,%s,%s,%s,%s',class_method,class_type,model,balance_triggers,balance_classes_method,bintrain,bintest,tfr_method,savelabels,labelsonly);
while strfind(str_settings,',,') str_settings = regexprep(str_settings,',,',','); end % remove duplicating commas
if str_settings(end) == ',' str_settings(end) = []; end % and remove trailing comma if present
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
    frequencies = '2:2:30';
end
if ischar(channels) && strcmpi(channels,'all')
    channels = 'ALL';
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