function FT_EEG = select_time_from_FT_EEG(FT_EEG,timebool)
% select time from FT_EEG struct
% moving away from the complicated fieldtrip implementations to keep things simple
% J.J.Fahrenfort, 2017

% trialdim = find(strcmpi(regexp(FT_EEG.dimord,'_','split'),'rpt')); -> use this if you want to implement reshaping
if ~strcmpi(FT_EEG.dimord,'rpt_chan_time')
    error('select_time_from_FT_EEG expects time as third dimension. Need to reshape, can easily be implemented if required');
end
FT_EEG.trial = FT_EEG.trial(:,:,timebool);
FT_EEG.time = FT_EEG.time(timebool);
