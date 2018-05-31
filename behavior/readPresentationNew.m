function log = readPresentationNew(fName,newNames,oldNames,joinEvents,eventDurationToNextEvent)
% function log = readPresentationNew(fName,newNames,oldNames,joinEvents,eventDurationToNextEvent)
% Reads the presentation logfile and recodes the values to a more suitable
% format. Picture events are recoded from oldNames to newNames.
% - newNames is a cell array of the new condition names. Conditions will be
%   coded using the index number of this array.
% - oldNames is an array containing the original condition numbers/names.
% - each cell in nameCodes can be a number array or a string. Codes
%   with exact matches to an element in the number array will be recoded to
%   the new value. Strings will be evaluated as regular expressions, where
%   any (partial) match will be recoded to the new value. See the function
%   regexp for more information on how to build regular expressions.
% - codes not specified in oldNames are removed
% - joinEvents allows one to specify whether trial codes of the (new) picture
%   events should be joined together when they are contiguously placed
%   in the same trial. Event durations are recoded to reflect this merge
%   operation.
% - if eventDurationToNextEvent is true, the duration of each picture event
%   in a trial is recoded to last until the next event (useful for blocked
%   designs).
%
% Example: 
% log = readLog('mylog.log',{'target', 'nontarget', 'mask'},{[1 2 3], [4 5],'[M m]ask\d'},true);
% recodes exact matches of 1, 2 or 3 to 'target', exact matches of 4 and 5
% to 'nontarget' and all matches containing the word 'Mask' (either with a
% capital M or a lower case m) followed by a digit between 0 and 9 to
% 'mask'. Contiguous trial events are joined together.

% By J.J.Fahrenfort, UvA 2012

if ~exist(fName)
    fName = getMatch(fName);
end
if numel(newNames)~=numel(oldNames)
    error('newNames and oldNames should be cell arrays of equal size');
end
if nargin < 4
    joinEvents = false;
end
if nargin < 5
    eventDurationToNextEvent = false;
end

% get basic info
fileInfo = dir(fName);
log.logFileCreated = fileInfo.date;
log.logFileName = fileInfo.name;
log.nameOfExperiment = fileInfo.name(strfind(fileInfo.name,'-')+1:regexp(fileInfo.name,'\d?\.')-1);
fid = fopen(fName,'r');
log.header1 = fgetl(fid);
log.header2 = fgetl(fid);
fgetl(fid);

% extract variable names
varNames = getEntries(fgetl(fid));
fgetl(fid);

% extract all values
cLine = 1;
dataline = fgetl(fid);
while dataline ~= -1
    entries = getEntries(dataline);
    nToExtract = min([numel(entries) numel(varNames)]);
    for c = 1:nToExtract
        data.(varNames{c}){cLine} = entries{c};
    end
    dataline = fgetl(fid);
    cLine = cLine + 1;
end

% recode values to doubles where possible and recode time to milliseconds
data.Trial = str2double(data.Trial);
timeColumns   = {'Time', 'TTime', 'Uncertainty', 'Duration', 'ReqTime', 'ReqDur'};
for c = 1:numel(timeColumns)
    data.(timeColumns{c}) = round(str2double(data.(timeColumns{c}))/10);
end

% get index array for each event type
log.subjectName = unique(data.Subject);
eventTypes = {'Pulse', 'Picture', 'Response', 'Quit', 'EyeEvent'};
for c = 1:numel(eventTypes)
    eventIndex.(eventTypes{c}) = find(strcmp(data.Event_Type,eventTypes{c}));
end

% shift time to either the first pulse or the first picture
if numel(eventIndex.Pulse > 0)
    startIndex = eventIndex.Pulse(1);
else
    startIndex = eventIndex.Picture(1);
end
data.Time = (data.Time - data.Time(startIndex));

% remove events prior to time 0
for c = 1:numel(eventTypes)
    eventIndex.(eventTypes{c})(eventIndex.(eventTypes{c})<startIndex) = [];
end

% create a new list of picture events, recoded from old events
% if old events are not defined, they are not included in the new list
event = [];
time = [];
duration = [];
oriCodes = [];
trial = [];
for cNew = 1:numel(newNames)
    criteria = oldNames{cNew};
    if isnumeric(criteria)
        tempIndex = [];
        criteria = strread(num2str(criteria),'%s');
        for c = 1:numel(criteria)
            tempIndex = [tempIndex find(strcmp(data.Code(eventIndex.Picture),criteria{c}))];
        end
    elseif ischar(criteria)
        tempIndex = find(~cellfun('isempty',regexp(data.Code(eventIndex.Picture),criteria)));
    else
        error('errr, your oldNames array contains content that cannot be evaluated');
    end
    event = [event repmat(cNew,1,numel(tempIndex))];
    time = [time data.Time(eventIndex.Picture(tempIndex))];
    duration = [duration data.Duration(eventIndex.Picture(tempIndex))];
    oriCodes = [oriCodes data.Code(eventIndex.Picture(tempIndex))];
    trial = [trial data.Trial(eventIndex.Picture(tempIndex))];
end

% join contiguous picture events if they are contained in the same trial
if joinEvents
    indexToRemove = [];
    trialNrs = unique(trial);
    for cTrials = 1:numel(trialNrs)
        indexOfTrial = find(trial == trialNrs(cTrials));
        codesInTrial = event(indexOfTrial);
        clusters = findClusters(codesInTrial,indexOfTrial);
        for c = 1:numel(clusters)
            if numel(clusters{c}) > 1
                duration(clusters{c}(1)) = time(clusters{c}(end)) + duration(clusters{c}(end)) - time(clusters{c}(1));
                indexToRemove = [ indexToRemove clusters{c}(2:end)];
            end
        end
    end
    % remove all obsolete events
    event(indexToRemove) = [];
    time(indexToRemove) = [];
    duration(indexToRemove) = [];
    oriCodes(indexToRemove) = [];
    trial(indexToRemove) = [];
end

% properly sort timeline
[log.stimuli.time timeIndex] = sort(time);
log.stimuli.event = event(timeIndex);
log.stimuli.duration = duration(timeIndex);
log.stimuli.oriCodes = oriCodes(timeIndex);
log.stimuli.trial = trial(timeIndex);
log.stimuli.conditionNames = newNames;

% recode each picture event to have a duration lasting until the next picture event
if eventDurationToNextEvent
    log.stimuli.duration(1:end-1) = log.stimuli.time(2:end)-log.stimuli.time(1:end-1);
end

% get other relevant data
log.pulse = data.Time(eventIndex.Pulse);
log.eye = data.Time(eventIndex.EyeEvent);
log.response.time = data.Time(eventIndex.Response);
log.response.event = str2double(data.Code(eventIndex.Response));

fclose all;

% function that matches which file has that pattern occuring first in its name
function [logfile] = getMatch(fName)
dirNames = dir;
dirNames = {dirNames(:).name};
matches = regexp(dirNames,fName);
indexOfMatches = find(~cellfun('isempty',matches));
if ~isempty(indexOfMatches)
    matchPos = cell2mat(matches);
    [ dummy bestMatch ] = min(matchPos);
    logfile = dirNames{indexOfMatches(bestMatch)}
else
    error([fName ' is not (part of) an existing filename in this directory.']);
end

% function to replace spaces by underscores and split up the string into a
% cell array using the tab delimiter 
function [entries] = getEntries(logLine)
logLine(regexp(logLine,' ')) = '_';
entries = regexp(logLine,'\t','split');

% function to determine clusters containing the index numbers of contiguous same-value code elements
function [clusterIndex] = findClusters(codes,codeIndex)
cluster = cumsum([1 diff(codes)~=0]);
clusterVals = unique(cluster);
for c = 1:numel(clusterVals)
    clusterIndex{c} = codeIndex(cluster==clusterVals(c));
end

