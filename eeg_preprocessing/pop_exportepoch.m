function events = pop_exportepoch(EEG,continuous)
% POP_EXPORTEPOCH exports epoch event values from an EEGLAB EEG structure. It only exports event
% values that occur at time 0 of your epoch, and it exports them as integers. These event values
% correspond to the values that the ADAM toolbox will be using. It is the function that ADAM uses
% internally to read EEGLAB events from your EEGLAB file. You can modify and re-import these event
% values using: EEG = pop_importepoch(EEG,events,{'eventtype'},'typefield','eventtype');
% If continuous data are exported, the second column of events contains time stamps expressed in
% milliseconds (!)
%
% Use as:
%   events = pop_exportepoch(EEG);
%
% internal function of the ADAM toolbox, by J.J.Fahrenfort, UvA/VU, 2018
%
% See also ADAM_MVPA_FIRSTLEVEL, READ_RAW_DATA
if nargin < 2
    continuous = false;
end

if continuous % extracting all continuous events, with time stamps
    events = NaN(numel(EEG.event),2);
    for cEvents = 1:numel(EEG.event)
        event = EEG.event(cEvents).type;
        if ischar(event)
            event = string2double(event);
            if numel(event) > 1
                warning('Cannot convert event to a single numeric value, taking the first value');
                event = event(1);
            end
            if isempty(event) || isnan(event)
                event = NaN;
                warning('Event value is NaN when converted to a numeral');
            end
        end
        events(cEvents,1) = event;
        pointlatency = EEG.event(cEvents).latency; % expressed in samples (!!!!)
        % eventlatency = round(pointlatency/EEG.srate*1000);
        eventlatency = EEG.times(round(pointlatency)); % expressed in ms. :-)
        events(cEvents,2) = eventlatency;
    end
else % extracting only epoched events that occur at time zero of each epoch, without time stamp
    events = NaN(numel(EEG.epoch),1);
    for cEvents = 1:numel(EEG.epoch)
        eventlatency = EEG.epoch(cEvents).eventlatency;
        if iscell(eventlatency)
            eventlatency = cell2mat(eventlatency);
        end
        timezero = find(eventlatency == 0);
        if isempty(timezero)
            error('Cannot find an event code with event latency 0, which is the critical event this function attempts to extract.');
        end
        event = EEG.epoch(cEvents).eventtype;
        if iscell(event)
            event = event{timezero};
        end
        if ~ischar(event) && numel(event) > 1
            event = event(timezero);
        end
        if ischar(event)
            event = string2double(event);
            if numel(event) > 1
                warning('Cannot convert event to a single numeric value, taking the first value');
                event = event(1);
            end
            if isempty(event) || isnan(event)
                event = NaN;
                warning('Event value is NaN when converted to a numeral');
            end
        end
        events(cEvents,1) = event;
    end
end
