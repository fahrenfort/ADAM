function FT_EEG = fix_dimord(FT_EEG)
% function FT_EEG = fix_dimord(FT_EEG)
% fixes the order of dimensions in FT_EEG to channel x time x trial
dims = regexp(FT_EEG.dimord, '_', 'split');
chandim = find(strcmp(dims,'chan'));
timedim = find(strcmp(dims,'time'));
trialdim = find(strcmp(dims,'rpt'));
if isempty(chandim | timedim | trialdim ) || numel(dims) ~= 3
    error('incorrect dimensions: should have time, channel and trial. If you also have frequency information, please use classify_TFR_data.');
end
if ~strcmp(FT_EEG.dimord,'chan_time_rpt')
    FT_EEG.trial = squeeze(permute(FT_EEG.trial,[chandim timedim trialdim ]));
    FT_EEG.dimord = 'chan_time_rpt';
end