function FT_EEG = eeglab2ft(EEG,filepath)
% function FT_EEG = eeglab2ft(EEG,filepath)
% Wrapper function to import EEG lab data to a fieldtrip, including
% conditions
% EEG can either be an eeglab struct, or be the filename which contains the
% eeglab data. In the latter case, filepath must also be specified
%
% J.J.Fahrenfort, VU 2014, 2016

if ~isstruct(EEG)
    if isempty(filepath)
        error('if you want to read in the data, you need to specify the path');
    end
    % remove .set if present
    EEG(strfind(EEG,'.set'):end) = [];
    % read data from EEGlab
    EEG = pop_loadset('filename',[EEG '.set'],'filepath',filepath);
end
FT_EEG = eeglab2fieldtrip(EEG, 'preprocessing');
% also get events
for cEvents = 1:numel(EEG.epoch)
    if iscell(EEG.epoch(cEvents).eventlatency)
        index = cell2mat(EEG.epoch(cEvents).eventlatency) == 0;
    else
        index = EEG.epoch(cEvents).eventlatency == 0;
    end
    event = EEG.epoch(cEvents).eventtype(index);
    if iscell(event)
        event = event{1};
    end
    if ischar(event)
        event = string2double(event);
        if numel(event) > 1
            warning('cannot convert event to a single numeric value, taking the first value');
            event = event(1);
        end
        if isempty(event) || isnan(event)
            event = NaN;
            warning('event value is NaN when converted to a numeral');
        end
    end
    FT_EEG.trialinfo(cEvents,1) = event;
end
clear EEG;

