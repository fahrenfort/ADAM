function [matObj, freq, time, trialinfo, chanindex, chandim, timedim, trialdim, freqdim] = read_mat_file(filename,channels)
matObj = matfile(filename);
dimord = matObj.dimord;
dims = regexp(dimord, '_', 'split');
chandim = find(strcmp(dims,'chan'));
timedim = find(strcmp(dims,'time'));
trialdim = find(strcmp(dims,'rpt'));
freqdim = find(strcmp(dims,'freq'));
if isempty(chandim | timedim | trialdim | freqdim) || numel(dims) ~= 4
    error('incorrect dimensions: should have time, channel, trial, and frequency. Either first compute TFRs or use classify_RAW_FT_data.');
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
% get electrode indices
new_channels = intersect(label,channels);
if isempty(channels)
    error(['cannot find any of the electrodes specified in ' channels]);
end
if numel(new_channels) < numel(channels)
    warning(['could not find all of the electrodes specified in ' channels]);
end
chanindex = zeros(1,numel(new_channels));
for cEl = 1:numel(new_channels)
    chanindex(cEl) = find(strcmp(new_channels{cEl},label));
end