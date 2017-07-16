function adam_MVPA_firstlevel(cfg)
% This function runs a first level analysis based on the parameters
% specificied in the cfg struct. Unspecified parameters take on default
% values. The following parameters can be specified:
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
channels = 'all';           % uses all electrodes except the EOG electrodes, other options: 'OCCIP' etc, see select_channels.m, only works for 64 channel EEG (not MEG)
nfolds = 10;                % trains on 90% (9/10) and tests on 10% (1/10)
crossclass = 'no';          % only trains and tests on the same time points
resample = 'no';            % does not resample the data
erp_baseline = 'no';        % [-.1,0] baselines from -100 to 0 ms
class_method = 'accuracy';  % computes classification accuracy (other options, e.g.: dprime, hr-far ( make sure that the first class is 'signal' and the second the 'noise')
class_type = 'linear';      % classifier type, e.g. 'linear' or 'diaglinear'
model = 'BDM';              % performs decoding rather than a forward encoding model
raw_or_tfr = 'raw';         % performs the analysis on the raw data rather than the time-frequeny data
balance_triggers = 'yes';   % balances triggers within each stimulus class, so that each class contains an equal amount of trigger values (discarding leftover triggers)
bintrain = 'no';            % average across triggers within a class on the training side
bintest = 'no';             % average across triggers witin a class on the testing side

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
if strcmpi(balance_triggers,'yes')
    balance_triggers = 'balanced';
else
    balance_triggers = 'unbalanced';
end
if strcmpi(bintrain,'yes')
    bintrain = 'bintrain';
else
    bintrain = '';
end
if strcmpi(bintest,'yes')
    bintrain = 'bintest';
else
    bintrain = '';
end
str_settings = sprintf('%s,%s,%s,%s,%s,%s',class_method,class_type,model,balance_triggers,bintrain,bintest);
% other settings
if strcmpi(crossclass,'no')
    crossclass = '0';
else
    crossclass = '1';
end
if strcmpi(resample,'no')
    resample = '0';
elseif ~ischar(resample)
    resample = num2str(resample);
end
crossclass_resample = sprintf('%s,%s',crossclass,resample);
if strcmpi(erp_baseline,'no')
    erp_baseline = '0';
end
if ischar(channels) && strcmpi(channels,'all')
    channels = '1';
end

% run analysis
for cSubj = 1:numel(filenames)
    %try
        if strcmpi(raw_or_tfr,'raw')
            classify_RAW_eeglab_data(datadir,filenames{cSubj},outputdir,nfolds,channels,str_settings,crossclass_resample,erp_baseline,class_spec{:});
        elseif strcmpi(raw_or_tfr,'TFR')
            disp('classify_TFR_from_eeglab_data still needs to be implemented in this function, is like 5 mins work, but not done yet');
        end
    %catch ME
    %    disp([ME.message ', skipping subject ' filenames{cSubj}]);
    %end
end