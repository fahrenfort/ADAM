function FT_EEG = select_channels_from_FT_EEG(FT_EEG,channelnames)
% select channels from FT_EEG struct and averages across them
% moving away from the complicated fieldtrip implementations to keep things simple
% J.J.Fahrenfort, 2017

% trialdim = find(strcmpi(regexp(FT_EEG.dimord,'_','split'),'rpt')); -> use this if you want to implement reshaping
if ~strcmpi(FT_EEG.dimord,'rpt_chan_time')
    error('select_channels_from_FT_EEG expects channels as second dimension. Need to reshape, can easily be implemented if required');
end
chans2keep = find(ismember(FT_EEG.label,channelnames));
if isempty(chans2keep)
    error(['Cannot find electrodes: ' cellarray2csvstring(channelnames)]);
end
FT_EEG.label = cellarray2csvstring(channelnames);
FT_EEG.trial = mean(FT_EEG.trial(:,chans2keep,:),2);
