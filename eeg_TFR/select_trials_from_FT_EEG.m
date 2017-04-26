function FT_EEG = select_trials_from_FT_EEG(FT_EEG,trialindex)
% select trials from FT_EEG struct
% moving away from the complicated fieldtrip implementations to keep things simple
% J.J.Fahrenfort, 2016

% trialdim = find(strcmpi(regexp(FT_EEG.dimord,'_','split'),'rpt')); -> use this if you want to implement reshaping
if ~strcmpi(FT_EEG.dimord,'rpt_chan_time')
    error('select_trials_from_FT_EEG expects trials as first dimension. Need to reshape, can easily be implemented if required');
end
FT_EEG.trial = FT_EEG.trial(trialindex,:,:);
FT_EEG.trialinfo = FT_EEG.trialinfo(trialindex);
