%% read files into EEGLAB and save the result
function EEG = readbva_into_eeglab(filename)
% read in data, subtract EOG, replace trigger values as needed, look up or
% insert channel info

% refchannels
refchannel = 'Ear_R';

% EOG above / below
horchannels = {'HEOG_R' 'HEOG_L'};

% eeglab path
% eeglab_path = [getenv('HOME') filesep 'Documents/matlab_toolboxes/eeglab13_4_4b'];
eeglab nogui;

% add eeglab path
% if (~isdeployed)
%     addpath(eeglab_path);
% end

ext = 'vhdr'; % 'xhdr'

% read data
if ~strcmpi(filename(end-3:end),ext)
    filename = [filename '.' ext];
end
[fpath, fname, ext ] = fileparts(filename);
if exist(filename,'file')
    % EEG = pop_fileio(filename);
    disp(filename);
    [EEG, com] = pop_loadbv(fpath, [fname ext]);
    %EEG = pop_loadbv(filename);
else
    error([fpath filesep fname ' does not seem to exist']);
end

% re-code triggers
for c=1:numel(EEG.event)
    code = EEG.event(c).type;
    code(regexp(code,'[S ]')) = [];
    EEG.event(c).type = code;
end

% deblank
for c=1:numel(EEG.chanlocs)
    EEG.chanlocs(c).labels = strtrim(EEG.chanlocs(c).labels);
end

% determine reference channel
ref = find(strcmpi(refchannel,{EEG.chanlocs.labels}));
if numel(ref) ~= 1
    error('could not find specified reference channel');
end

% re-reference to leftover earlobe
EEG.data = EEG.data - repmat(EEG.data(ref, :) * 0.5, [EEG.nbchan 1]);
EEG.data = EEG.data(setdiff(1:EEG.nbchan,ref),:);
EEG.chanlocs = EEG.chanlocs(setdiff(1:EEG.nbchan,ref));
EEG.nbchan = EEG.nbchan - 1;

% % fix triggers if needed
% % EEG.event(x).type = condition
% % EEG.event(x).latency = sample
% if replace_triggers
%     EEG.event = EEG.event(1:numel(triggers.(filename)));
%     [EEG.event.type] = deal(triggers.(filename).value);
%     [EEG.event.latency] = deal(triggers.(filename).sample);
% end

% % re-reference ocular channels
% hor = [ find(strcmpi(horchannels{1},{EEG.chanlocs.labels})) find(strcmpi(horchannels{2},{EEG.chanlocs.labels}))];
% if numel(hor) ~= 2 
%     warning('cannot find horizontal eye channels');
% else
%     horEyeData = EEG.data(hor(1),:)-EEG.data(hor(2),:);
%     % remove one
%     EEG = pop_select( EEG,'nochannel',horchannels(1));
%     % and replace the other with the substraction
%     hor2indx = find(strcmpi(horchannels{2},{EEG.chanlocs.labels}));
%     EEG.data(hor2indx,:) = horEyeData;
%     EEG.chanlocs(hor2indx).labels = 'HEOG';
% end

% some bookkeeping
if isfield(EEG.chaninfo,'nodatchans')
    EEG.chaninfo = rmfield(EEG.chaninfo,'nodatchans');
end

EEG = eeg_checkset(EEG);
