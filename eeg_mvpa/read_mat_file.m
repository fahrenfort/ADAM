function [matObj, freq, time, trialinfo, chanindex, chandim, timedim, trialdim, freqdim] = read_mat_file(filename,channels)
matObj = matfile(filename);
dimord = matObj.dimord;
dims = regexp(dimord, '_', 'split');
chandim = find(strcmp(dims,'chan'));
timedim = find(strcmp(dims,'time'));
trialdim = find(strcmp(dims,'rpt'));
freqdim = find(strcmp(dims,'freq'));
if isempty(chandim | timedim | trialdim | freqdim) || numel(dims) ~= 4
    error('incorrect dimensions: should have time, channel, trial, and frequency.');
end
freq = matObj.freq;
label = matObj.label;
time = matObj.time;
vars = whos('-file',filename);
if ismember('trialinfo', {vars.name})
    trialinfo = matObj.trialinfo;
else
    trialinfo = [];
end
% get electrode indices: chanindex returns the indices for which the order of the channels in
% matObj.data(chanindex) is the same as in the input parameter 'channels'
[new_channels,~,chanindex] = intersect(channels,label,'stable'); 
if numel(new_channels) ~= numel(channels)
    disp(['WARNING: could not find all of the electrodes specified in ' cellarray2csvstring(channels)]);
end