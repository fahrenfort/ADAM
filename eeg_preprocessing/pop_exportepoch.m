function events = pop_exportepoch(EEG)
% exports epoch triggers that occur at time 0, as integers
events = zeros(1,numel(EEG.epoch));
for c=1:numel(EEG.epoch);
    eventlatency = EEG.epoch(c).eventlatency;
    if iscell(eventlatency)
        eventlatency = cell2mat(eventlatency);
    end
    eventtype = EEG.epoch(c).eventtype;
    if iscell(eventtype)
        event = eventtype{eventlatency==0};
    else
        event = eventtype(eventlatency==0);
    end
    if ischar(event)
        event = string2double(event);
    end
    events(c) = event;
end
