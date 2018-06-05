function out = addResp2Trials(data, time_window, eventvalues, last)
% addResp2Trials adds responses to the stimulus field after reading in presentation logfiles using
% readPresentationNew.
%
% inputs:
%           time_window  -  within which response has to take place in ms. (default: 200-1500)
%           eventvalues  -  the response event values that are included in the operation (default:
%                           empty, reading in all responses)
%           last        -   if true, takes the last response, otherwise takes the first response
%                           (default: false)
% outputs:
%           out has the same structure as data, but now with response events and their associated
%           reaction times added to the stimulus field
%
% By J.J.Fahrenfort, UvA/VU 2018

if nargin < 4
    last = false;
end
if nargin < 3
    eventvalues = [];
end
if nargin < 2
    time_window = [200, 1500];
end

% add responses to stimulus table
out = data;
for cStim = 1:numel(data.stimuli.time)
    startwin = data.stimuli.time(cStim)+time_window(1);
    stopwin = data.stimuli.time(cStim)+time_window(2);
    if isempty(eventvalues)
        respIndex = find(data.response.time >= startwin & data.response.time <= stopwin);
    else
        respIndex = find(data.response.time >= startwin & data.response.time <= stopwin & ismember(data.response.event,eventvalues));
    end
    if ~isempty(respIndex)
        responses = data.response.event(respIndex);
        RT = data.response.time(respIndex) - data.stimuli.time(cStim);
    else
        responses = NaN;
        RT = NaN;
    end
    if last
        out.stimuli.response(cStim) = responses(end);
        out.stimuli.RT(cStim) = RT(end);
    else
        out.stimuli.response(cStim) = responses(1);
        out.stimuli.RT(cStim) = RT(1);      
    end
end