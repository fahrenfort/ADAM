function [matObj, dim_params] = read_mat_file(filename,channels,frequency,timeindex)
% creates pointer to FT_EEG data and computes a dim_params.index definition that allows one to
% access only a part of the data, as defined in channels, frequency and/or timeindex.

if nargin<4
    timeindex = [];
end
if nargin<3
    frequency = [];
end
if nargin<2
    channels = [];
end

% this is the file
matObj = matfile(filename);

% these are the dimensions in the file
dimord = matObj.dimord;
dims = regexp(dimord, '_', 'split');
chandim = find(strcmp(dims,'chan'));
timedim = find(strcmp(dims,'time'));
trialdim = find(strcmp(dims,'rpt'));
freqdim = find(strcmp(dims,'freq'));

% make cell indexing to extract data
index = cell(1,numel(dims));
index(:) = {':'};

% frequency
if ~isempty(freqdim)
    freq = matObj.freq;
    actualfrequencies = round(freq*100)/100;
    % get frequency index
    if ~isempty(frequency)
        freqindex = find(frequency==actualfrequencies);
        if isempty(freqindex)
            disp('WARNING: cannot find an exact match for that frequency, attempting to find the closest match');
            freqindex = nearest(actualfrequencies,frequency);
            if isempty(freqindex)
                error(['error, cannot find a matching frequency for frequency ' num2str(frequency) ', giving up now']);
            end
        end
        index{freqdim} = freqindex;
    end
end

% channels
if ~isempty(chandim)
    label = matObj.label;
    % get electrode indices in the same order as in the input parameter 'channels'
    if ~isempty(channels)
        [new_channels,~,chanindex] = intersect(channels,label,'stable');
        if numel(new_channels) ~= numel(channels)
            disp(['WARNING: could not find all of the electrodes specified in ' cellarray2csvstring(channels)]);
        end
        index{chandim} = chanindex;
    end
end

% time
if ~isempty(timedim) && ~isempty(timeindex)
    index{timedim} = timeindex;
end

% return relevant info
dim_params.chandim = chandim;
dim_params.timedim = timedim;
dim_params.trialdim = trialdim;
dim_params.freqdim = freqdim;
dim_params.index = index;
