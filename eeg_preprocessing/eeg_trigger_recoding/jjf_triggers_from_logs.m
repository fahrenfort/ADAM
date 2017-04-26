function [ output, log ] = jjf_triggers_from_logs(datafile,logfile)
% function to read event triggers and place them in event struct
% By JJ Fahrenfort, 2014

% READ event and header from EEG
event = ft_read_event(datafile);
header = ft_read_header(datafile);

% RECOVER STIM values (binary values from status channel)
valueCell = {event.value};
emptyInd = find(cellfun(@isempty,valueCell));
for c = 1:numel(emptyInd)
    event(emptyInd(c)).value = 999; % replace empty values with 999
end
binVal = dec2bin([event.value]);
numBits = numel(binVal(1,:));
% get relevant values
triggers = bin2dec(binVal(:,(numBits-6):(numBits-0)));
samples = [event.sample];
sampletimes = jjf_sample2time(samples,header.Fs); 

% REMOVE values prior to CM_in_range
CM_index = find(strcmp('CM_in_range',{event.type}));
indexOfFirstStim = CM_index + 1;

% make some corrections due to acquisition error
if sum(strfind(datafile,'WM_pilot_pp1s2.bdf'))
    fprintf('making corrections! : started recording too late WM_pilot_pp1s2.bdf\n');
    indexOfFirstStim = find(triggers==72,1); % 72 was the first trigger code that was sent in a trial
    usefultriggers = numel(triggers(indexOfFirstStim:end));
elseif sum(strfind(datafile,'WM_pilot_pp1s2_2.bdf'))
    fprintf('making corrections! this session was broken off: WM_pilot_pp1s2_2.bdf\n');
    indexOfFirstStim = find(triggers==72,1); % 72 was the first trigger code that was sent in a trial
    usefultriggers = numel(triggers(indexOfFirstStim:end));
elseif sum(strfind(datafile,'WM_pilot_pp2s1.bdf'))
    fprintf('making corrections! : started recording too late WM_pilot_pp2s1.bdf\n');
    indexOfFirstStim = find(triggers==72,1); % 72 was the first trigger code that was sent in a trial
    usefultriggers = numel(triggers(indexOfFirstStim:end));
elseif sum(strfind(datafile,'WM_pilot_pp4s1.bdf'))
    fprintf('making corrections! there is one trigger too many in this dataset: WM_pilot_pp4s1.bdf\n')
    indexOfWrong = find(triggers==20);
    indexOfWrong = indexOfWrong(end);
    triggers(indexOfWrong) = [];
    samples(indexOfWrong) = [];
    sampletimes(indexOfWrong) = [];
elseif sum(strfind(datafile,'WM_pilot_pp4s2.bdf'))
    fprintf('making corrections! : started recording too late WM_pilot_pp4s2.bdf\n');
    indexOfFirstStim = find(triggers==72,1); % 72 was the first trigger code that was sent in a trial
    usefultriggers = numel(triggers(indexOfFirstStim:end));
end
% and continue as normal

% EXTRACT events and timing from BDF
origtriggers = triggers(indexOfFirstStim:end);
stimsamples = samples(indexOfFirstStim:end);
stimsampletimes = sampletimes(indexOfFirstStim:end);
nTriggers = numel(origtriggers);

% EXTRACT events from logfile
logtable = readtable(logfile);
stimIndex = find(~isnan(logtable.condition));
nTriggersInLog = numel(stimIndex);
% make some more corrections if necessary 
if exist('usefultriggers')
    if sum(strfind(datafile,'WM_pilot_pp1s2.bdf'))
        stimIndex = stimIndex(3:usefultriggers/4+3); % here first part of condition log is taken
    else
        stimIndex = stimIndex(end-(usefultriggers/4)+1:end); % here last part of condition log is taken
    end
end
log.condition = logtable.condition(stimIndex);
log.tasktype = logtable.tasktype(stimIndex);
log.setsize = logtable.setsize(stimIndex);
log.block_nr = logtable.block_nr(stimIndex);
log.response = logtable.key_resp_2_corr(stimIndex);
log.resp_rt = round(1000*(logtable.key_resp_2_rt(stimIndex)));% - logtable.search_target_started(stimIndex)));
log.bar_duration = round(1000*(logtable.retention_started(stimIndex) - logtable.bar_started(stimIndex)));
log.retention_duration = round(1000*(logtable.search_target_started(stimIndex) - logtable.retention_started(stimIndex)));


% STIMULUS CONDITION RECODING AS WAS DONE IN BIOSEMI (USED FOR CHECKING)
% tasktype 1 = orientation (add 0), tasktype 2 = color (add 100 to condition)
condition = log.condition;
condition(find(log.tasktype==2)) = condition(find(log.tasktype==2)) + 100;
condition(find(log.setsize==10)) = condition(find(log.setsize==10)) + 50;
% BUT ACTUALLY, SET SIZE DID NOT WORK IN BEHAVIOR, so let's not include it
% and stick to coding for color and orientation
cond = log.condition;
cond(find(log.tasktype==2)) = cond(find(log.tasktype==2)) + 100;

% INTEGRATE TRIGGERS from logfile into EEG
% there are 4 trigger moments: 
% (1) onset bar (add 1000 to condition)
% (2) onset retention interval (add 0 to condition)
% (3) onset memory array (add 2000 to condition)
% (4) response (add 3000 to condition)
indexOfBar = 1:4:numel(cond)*4;
indexOfRet = 2:4:numel(cond)*4;
indexOfMem = 3:4:numel(cond)*4;
indexOfResp = 4:4:numel(cond)*4;
oldstimtriggers(indexOfBar) = 200;
oldstimtriggers(indexOfRet) = condition;
oldstimtriggers(indexOfMem) = 250;
oldstimtriggers(indexOfResp) = log.response+254;
newstimtriggers(indexOfBar) = cond+1000;
newstimtriggers(indexOfRet) = cond;
newstimtriggers(indexOfMem) = cond+2000;
newstimtriggers(indexOfResp) = log.response+3000;
% and add in values for other noteworthy information
block_nrs(indexOfBar) = log.block_nr;
block_nrs(indexOfRet) = log.block_nr;
block_nrs(indexOfMem) = log.block_nr;
block_nrs(indexOfResp) = log.block_nr;
setsizes(indexOfBar) = log.setsize;
setsizes(indexOfRet) = log.setsize;
setsizes(indexOfMem) = log.setsize;
setsizes(indexOfResp) = log.setsize;
responses(indexOfBar) = log.response;
responses(indexOfRet) = log.response;
responses(indexOfMem) = log.response;
responses(indexOfResp) = log.response;

% CHECK AND WARN IF number of triggers in EEG is too low or too high
if nTriggers < nTriggersInLog*4
    fprintf([ '!! number of triggers in EEG (' num2str(nTriggers) ') is smaller than number of triggers in logfile (' num2str(nTriggersInLog*4) '):\n' datafile '\n']);
elseif nTriggers > nTriggersInLog*4
    error('number of triggers in EEG larger than number of triggers in the logfile -> no use, give up');
end

% FINAL CHECK TO SEE IF ORIGINAL VALUES MATCH
[correct, outcome] = jjf_checktriggers(origtriggers(2:4:end),condition);
if correct
    fprintf(['Oh happy joy, all ' num2str(numel(outcome)) ' trigger values match!\n\n']);
else
    fprintf('not looking good: ');
    fprintf([num2str(sum(~outcome)) ' triggers of ' num2str(numel(outcome)) ' trigger values do not match :-(\n\n']);
end

% COPY THE VALUES BACK TO OUTPUT
% (we've thrown out the CM_in_range and triggers prior to that)
for c = 1:nTriggers
    if c>=2 && ~jjf_checktriggers(origtriggers(c),oldstimtriggers(c))
        error('somethings wrong after all');
    end
    output(c).sample = stimsamples(c);
    output(c).value = newstimtriggers(c);
    output(c).origvalue = origtriggers(c);
    output(c).type = 'STATUS';
    output(c).offset = stimsampletimes(c);
    output(c).block_nr = block_nrs(c);
    output(c).setsize = setsizes(c);
    output(c).response = responses(c);
end

