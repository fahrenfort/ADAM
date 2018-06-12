function adam_create_group_training_sets(cfg)
% adam_create_group_training_sets creates a 'group average' training set for each subject, which
% replaces each trial with a group average across all subjects, excluding the subject for which the
% training set is intended. Outputs these training sets in FT_EEG format.
%
% Inputs:
%   datadir         - string specifying the path to the input files (eeglab .set/.fdt or fieldtrip
%                     .mat), required.
%   filenames       - cell array of filenames for which the group training sets are created,
%                     required.
%   erp_baseline    - [int int] a baseline period specified in seconds, e.g. [-.1 0] or
%                     'no' (not recommended), required.
%   outdir          - where the results are saved, optional. If outdir is not specified, files are
%                     saved in datadir.
%   eventvals       - eventvals for which the operation is performed (if left empty, it will use all
%                     eventvals on which the data was epoched, optional.
%   keepratio       - .75 (default), the minimum number of subjects on which each trial average is
%                     based. If the trial average is based on fewer subjects than specified by this
%                     ratio, the trial will not be included in the training set, optional.
%   skipuntil       - 'subjectname' (optional field, string), if the function previously crashed and
%                     you want to skip until this subject, optional.
%
% A number of additional parameters can be specified, see read_raw_data for parameter values.
%
% Outputs: data is saved to datadir or outdir. Each filename is prepended with the prefix
% 'grptrain_'
%
% J.J.Fahrenfort, UvA/VU 2018

% for this function
datadir = [];
filenames = [];
outputdir = [];
erp_baseline = [];
skipuntil = [];

% unpack settings
v2struct(cfg);
missingparams = cellfun(@isempty,{datadir, filenames, erp_baseline});
if any(missingparams)
    paramstrings = {'datadir', 'filenames', 'erp_baseline'};
    error(['Missing cfg fields: ' cell2csv(paramstrings(missingparams),true)]);
end

% check existence of input dir
if ~exist(datadir,'dir')
    error([datadir ' does not exist.']);
end

% create target directory
if isempty(outputdir)
    outputdir = datadir;
elseif ~exist(outputdir,'dir')
    mkdir(outputdir);
end

% determine first subject
if isempty(skipuntil)
    startsubject = 1;
else
    startsubject = find(strcmpi(filenames,skipuntil));
    disp(['starting with subject ' filenames{startsubject}]);
end

for cFiles=startsubject:numel(filenames)
    cfg.files2use = setdiff(filenames,filenames{cFiles});
    FT_EEG = create_training_set(cfg);
    save(fullfile(outputdir,['grptrain_' filenames{cFiles}]),'FT_EEG','-v7.3');
end

function FT_EEG_OUT = create_training_set(cfg)
% settings for read data:
detrend_eeg = false;
shuffle_trials = true;
clean_window = [];
clean_data = false;
do_csd = false;
resample_eeg = false;
channelset = 'all';

% some more defaults
eventvals = [];
keepratio = .75;

% unpack
v2struct(cfg);
msettings = v2struct(detrend_eeg,shuffle_trials,clean_window,clean_data,do_csd,resample_eeg,channelset,erp_baseline);

% which files to use
filenames = files2use;

% read in first subject as basis
FT_EEG_ALL = read_raw_data(datadir,filenames{1},outputdir,msettings);

% determine size
time = FT_EEG_ALL.time;
label = FT_EEG_ALL.label;
nChans =  numel(label);
nTimepoints =  numel(time);

% which events to include?
if isempty(eventvals)
    eventvalsbase = unique(FT_EEG_ALL.trialinfo);
else
    eventvalsbase = eventvals;
end

% create event bins
for cEv = 1:numel(eventvalsbase)
    eventbin{cEv} = FT_EEG_ALL.trial(ismember(FT_EEG_ALL.trialinfo,eventvalsbase(cEv)),:,:);
    trialCounter{cEv} = ones(size(eventbin{cEv},1),1);
end

% keep memory low
clear FT_EEG_ALL;

% keep on adding the other subjects
for cFiles=2:numel(filenames)
    disp(['loading subject: ' filenames{cFiles}]);
    FT_EEG = read_raw_data(datadir,filenames{cFiles},outputdir,cfg);
    
    % check if sizes are the same
    if any([nChans~=size(FT_EEG.trial,2) nTimepoints~=size(FT_EEG.trial,3)])
        error(['not the same number of channels and/or samples in one of the data sets:' filenames{cFiles}]);
    end
    
    % check if channel ordering is the same
    if ~all(time == FT_EEG.time) || ~all(strcmp(label, FT_EEG.label))
        error(['Time points or channel (ordering) not the same in data:' filenames{cFiles}])
    end
    
    % which events to include?
    if isempty(eventvals)
        eventvals2add = unique(FT_EEG.trialinfo);
    else
        eventvals2add = eventvals;
    end
    
    % make list of all eventvals, per dataset
    for cEv = 1:numel(eventvals2add)
        eventbin2add = FT_EEG.trial(ismember(FT_EEG.trialinfo,eventvals2add(cEv)),:,:);
        baseIndex = find(eventvalsbase == eventvals2add(cEv));
        % if event was not present in base set, add it
        if isempty(baseIndex)
            N2add = size(eventbin2add,1);
            eventvalsbase(end+1) = eventvals2add(cEv);
            eventbin{end+1} = eventbin2add;
            trialCounter{end+1} = ones(N2add,1);
            % if more events exist in the new eventbin2add, increase size of eventbin
        else
            N2add = size(eventbin2add,1);
            Nbase = size(eventbin{baseIndex},1);
            if N2add > Nbase
                % increase eventbin
                eventbin{baseIndex} = vertcat(eventbin{baseIndex}, zeros(N2add-Nbase, nChans, nTimepoints));
                trialCounter{baseIndex} = vertcat(trialCounter{baseIndex}, zeros(N2add-Nbase,1));
                trialCounter{baseIndex} =  trialCounter{baseIndex} + 1;
            elseif N2add <= Nbase
                % increase eventbin2add
                eventbin2add = vertcat(eventbin2add, zeros(Nbase-N2add, nChans, nTimepoints));
                trialCounter{baseIndex}(1:N2add) =  trialCounter{baseIndex}(1:N2add) + 1;
            end
            eventbin{baseIndex} = eventbin{baseIndex} + eventbin2add;
        end
    end
end

% keep memory low
clear FT_EEG;

% create trialinfo for each bin
for cEv = 1:numel(eventvalsbase)
    trialinfo{cEv} = ones(size(eventbin{cEv},1),1)*eventvalsbase(cEv);
end

% merge eventbins into single trial and trialinfo and counter variables
trialCounter = vertcat(trialCounter{:});
trial = vertcat(eventbin{:});
trialinfo = vertcat(trialinfo{:});

% average!
trial = trial./repmat(trialCounter,[1 nChans nTimepoints]);

% prune the final trials based on keepratio
realratio = trialCounter/numel(filenames);
index2remove = realratio < keepratio;
trial(index2remove,:,:) = [];
trialinfo(index2remove) = [];

% return results
FT_EEG_OUT.time = time;
FT_EEG_OUT.label = label;
FT_EEG_OUT.trial = trial;
FT_EEG_OUT.trialinfo = trialinfo;
FT_EEG_OUT.dimord = 'rpt_chan_time';
