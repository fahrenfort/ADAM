function EEG = add_channel_info(EEG)
% double checks whether channel information is present, add if not
% J.J.Fahrenfort, VU 2017

% double check whether channel location information is present
channelnames = {EEG.chanlocs(:).labels};
nopos_channels = [];
for cEl=1:length(EEG.chanlocs)
    if (any(isempty(EEG.chanlocs(1,cEl).X)&isempty(EEG.chanlocs(1,cEl).Y)&isempty(EEG.chanlocs(1,cEl).Z)&isempty(EEG.chanlocs(1,cEl).theta)&isempty(EEG.chanlocs(1,cEl).radius)))
        nopos_channels = [nopos_channels cEl];
    end
end

% find EEG channels
eeg_channels = select_channels(channelnames,'EEG');
eeg_eog_channels = select_channels(channelnames,'EEG_EOG');
if isempty(eeg_channels)
    error('stopping now, there are no EEG channels in this set??');
end

% in case some EEG channels do not have location information
if any(ismember(eeg_channels,nopos_channels))
    disp(['WARNING: Channels ' num2str(nopos_channels) ' have incomplete location information. Now attempting to read in location information.']);
    try 
        EEG = pop_chanedit(EEG, 'lookup', findcapfile);
    catch
        try
            EEG = pop_chanedit(EEG, 'lookup', which('standard-10-5-cap385.elp'));
        catch
            error('Cannot retrieve location information for some channels. Fix this first!');
        end
    end
end