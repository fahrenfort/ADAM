function FT_EEG = detrend_FT_EEG(FT_EEG)
% detrends FT_EEG data
% by J.J.Fahrenfort, VU 2017
disp('Detrending FT_EEG data');
cfg = [];
FT_EEG = ft_timelockbaseline(cfg,FT_EEG);
dims = regexp(FT_EEG.dimord, '_', 'split');
chandim = find(strcmp(dims,'chan'));
timedim = find(strcmp(dims,'time'));
trialdim = find(strcmp(dims,'rpt'));
origsize = size(FT_EEG.trial);
FT_EEG.trial = reshape(permute(FT_EEG.trial, [timedim chandim trialdim]), [size(FT_EEG.trial,timedim) size(FT_EEG.trial,chandim)*size(FT_EEG.trial,trialdim)]);
FT_EEG.trial = detrend(FT_EEG.trial);
FT_EEG.trial = permute(reshape(FT_EEG.trial, [origsize(timedim) origsize(chandim) origsize(trialdim)]), [3 2 1]);
FT_EEG.dimord = 'rpt_chan_time';

% EEG.data = reshape(permute(EEG.data, [2 1 3]), [EEG.pnts EEG.nbchan * EEG.trials]);
% EEG.data = detrend(EEG.data);
% EEG.data = permute(reshape(EEG.data, [EEG.pnts EEG.nbchan EEG.trials]), [2 1 3]);