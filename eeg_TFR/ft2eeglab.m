function EEG = ft2eeglab(FT_EEG)
% This function converts an FT_EEG structure (the Fieldtrip timelock datatype, as used inside the
% ADAM toolbox) to an EEGLAB structure that can be read in and used inside EEGLAB.
% Internal function of ADAM toolbox, J.J.Fahrenfort, UvA/VU 2019

if ~isstruct(FT_EEG)
    error('This function takes a Fieldtrip struct as input.');
end

FT_EEG = fix_dimord(FT_EEG,'chan_time_rpt');

for c=1:numel(FT_EEG.label)
    chanlocs(c).labels = FT_EEG.label{c};
end
EEG.chanlocs   = pop_chanedit(chanlocs, 'lookup', trycapfile);
EEG.data       = single(FT_EEG.trial); % electrodes x time x trials
EEG.setname    = 'from FT_EEG';
EEG.filename   = '';
EEG.filepath   = '';
EEG.subject    = '';
EEG.group      = '';
EEG.condition  = '';
EEG.session    = [];
EEG.comments   = 'preprocessed with fieldtrip and/or ADAM';
EEG.nbchan     = size(FT_EEG.trial,1);
EEG.trials     = size(FT_EEG.trial,3);
EEG.pnts       = size(FT_EEG.trial,2);
EEG.srate      = FT_EEG.fsample;
EEG.xmin       = FT_EEG.time(1)/1000;    % assumes time is in seconds
EEG.xmax       = FT_EEG.time(end)/1000;
EEG.times      = FT_EEG.time; % in miliseconds
EEG.ref        = []; %'common';
EEG.event      = [];
EEG.epoch      = [];
EEG.icawinv    = [];
EEG.icasphere  = [];
EEG.icaweights = [];
EEG.icaact     = [];
EEG.saved      = 'no';
% also import event structure
EEG = pop_importepoch(EEG,FT_EEG.trialinfo(:,1),{'eventtype'},'typefield','eventtype');
EEG = eeg_checkset(EEG);