%% read files into EEGLAB and save the result
function EEG = readbdf_into_eeglab(filename)
% read in data, subtract EOG, replace trigger values as needed, look up or
% insert channel info

% some specifics
refchannels = {'EXG5' 'EXG6'};
% EOG above / below
verchannels = {'EXG3' 'EXG4'}; % verchannels = {'VEOGsT' 'VEOGsT'}; 
% EOG left / right
horchannels = {'EXG2' 'EXG1' }; % horchannels = {'HEOGsL' 'HEOGsR'}; 
% eeglab path
eeglab_path = [getenv('HOME') filesep 'Documents/matlab_toolboxes/eeglab13_4_4b'];

% add eeglab path
if (~isdeployed)
    addpath(eeglab_path);
end

% load triggers
datapath = filename(1:max(strfind(filename,filesep))-1);
if exist([datapath filesep 'triggercodes.mat'],'file');
    replace_triggers = 1;
    load([datapath filesep 'triggercodes'], 'triggers');
    disp('triggers will be replaced!');
else
    replace_triggers = 0;
end

% read data
if ~strcmpi(filename(end-3:end),'.bdf')
    filename = [filename '.bdf'];
end
if exist(filename,'file')
    disp(filename);
    %EEG = pop_biosig(filename);
    EEG = pop_readbdf_triggerfix(filename,[],73); % 73/80 refers to the status channel 
    % the tailor made pop_readbdf_triggerfix can be found in my /matlab_scripts/eeg_preprocessing folder
else
    error([datapath filesep filename ' does not seem to exist']);
end

% fix labels if needed
if exist([datapath filesep '68ChanLocsBiosemi.mat'],'file')
    load([datapath filesep '68ChanLocsBiosemi.mat']);
    load([datapath filesep '68ChanInfoBiosemi.mat']);
    EEG.chanlocs = chanlocs;
    EEG.chaninfo = chaninfo;
else
    disp('not reading external channel file, are you sure the labels are correct?');
end

% determine reference channels
ref = [find(strcmpi(refchannels{1},{EEG.chanlocs.labels})) find(strcmpi(refchannels{2},{EEG.chanlocs.labels}))];
if numel(ref) ~= 2
    error('could not find specified reference channels');
end

% re-reference
EEG = pop_reref(EEG, ref, 'refstate','averef');

% fix triggers if needed
% EEG.event(x).type = condition
% EEG.event(x).latency = sample
if replace_triggers
    EEG.event = EEG.event(1:numel(triggers.(filename)));
    [EEG.event.type] = deal(triggers.(filename).value);
    [EEG.event.latency] = deal(triggers.(filename).sample);
end

% re-reference ocular channels
ver = [ find(strcmpi(verchannels{1},{EEG.chanlocs.labels})) find(strcmpi(verchannels{2},{EEG.chanlocs.labels}))];
hor = [ find(strcmpi(horchannels{1},{EEG.chanlocs.labels})) find(strcmpi(horchannels{2},{EEG.chanlocs.labels}))];
if numel(hor) ~= 2 || numel(ver) ~= 2
    warning('cannot find horizontal and/or vertical eye channels');
else
    verEyeData = EEG.data(ver(1),:)-EEG.data(ver(2),:);
    horEyeData = EEG.data(hor(1),:)-EEG.data(hor(2),:);
    EEG.data(ref(1),:) = verEyeData;
    EEG.data(ref(2),:) = horEyeData;
    EEG.chanlocs(ref(1)).labels = 'VEOG';
    EEG.chanlocs(ref(2)).labels = 'HEOG';
    EEG.chanlocs(ver(1)).labels = 'SO1';
    EEG.chanlocs(ver(2)).labels = 'IO1';
    EEG.chanlocs(hor(1)).labels = 'LO1';
    EEG.chanlocs(hor(2)).labels = 'LO2';
    % remove unwanted channels
    EEG.data = EEG.data(1:max([ref ver hor]),:);
    EEG.chanlocs = EEG.chanlocs(1:max([ref ver hor]));
    EEG.nbchan = max([ref ver hor]);
end

% some bookkeeping
if isfield(EEG.chaninfo,'nodatchans')
    EEG.chaninfo = rmfield(EEG.chaninfo,'nodatchans');
end
EEG = pop_select( EEG,'nochannel',{ 'EXG1' 'EXG2' 'EXG3' 'EXG4' 'EXG5' 'EXG6' 'EXG7' 'EXG8' });
EEG = pop_chanedit(EEG, 'lookup',[eeglab_path filesep 'plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp']);
EEG = pop_editset(EEG, 'setname', filename(max(strfind(filename,filesep))+1:end-4));