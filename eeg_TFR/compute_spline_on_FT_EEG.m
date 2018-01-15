function FT_EEG_SPLINE = compute_spline_on_FT_EEG(FT_EEG,relevant_window_in_ms,field,basefreq)
% function FT_EEG_SPLINE = compute_spline_on_FT_EEG(FT_EEG,relevant_window_in_ms,field,basefreq)
% turn all ERPs into splines
% use 8hz as a base frequency
% relevant_window_in_ms sets the window to detect the peak of the spline in
% milliseconds (defaults to [50,750]).
% Johannes Fahrenfort, VU, 2016
npadding = round(FT_EEG.fsample/4); % 250 ms
if nargin < 4
    basefreq = 8;
end
if nargin < 3 || isempty(field)
    field = 'trial';
end
if nargin < 2
    relevant_window_in_ms = [50 750];
end
if ~isfield(FT_EEG,'dimord')
    error('the input dataset is not in the required standard fieldtrip format. you might want to run ft_timelockbaseline on FT_EEG to resolve this.');
else
    dims = regexp(FT_EEG.dimord, '_', 'split');
    chandim = find(strcmp(dims,'chan'));
    timedim = find(strcmp(dims,'time'));
    trialdim = find(strcmp(dims,'rpt'));
end
if isempty(chandim | timedim | trialdim) || numel(dims) > 3
    error('incorrect dimensions: should have time, channel and trial, cannot compute ERP on FT_EEG');
end
if ~isfield(FT_EEG,field)
    error([field ' should contain the data, but ' field ' is not a field of FT_EEG']);
end
indata = permute(FT_EEG.(field),[trialdim chandim timedim]);
outdata = zeros(size(indata));
window_indx = nearest(FT_EEG.time*1000,relevant_window_in_ms);

for cTrial = 1:size(indata,1)
    for cChan = 1:size(indata,2)
        data = squeeze(indata(cTrial,cChan,:));
        window_data = data(window_indx(1):window_indx(2));
        % get the max point
        if max(window_data) > abs(min(window_data)) % center around max or min value in the time array
            [~, peakindx] = max(window_data); 
        else
            [~, peakindx] = min(window_data); 
        end
        peakindx = peakindx(1); % pick the first one
        % continue
        paddata = padarray(data,npadding,'symmetric');
        peakindx = peakindx + npadding + window_indx(1) - 1; % new peakindex based on padded data
        % resample to max 8 hz
        newstepsize = FT_EEG.fsample/basefreq;
        % get new indexes based on this sampling frequency
        indx = round(1:newstepsize:numel(paddata));
        peakdif = indx(nearest(indx,peakindx))-peakindx; % center around max (or min) of trial to get the peak out
        newindx = indx - peakdif;
        newindx = newindx(newindx>0);
        newindx = newindx(newindx<=numel(paddata));
        if size(newindx) == 0
            error('this cannot happen');
        end
        spline_erp = spline(newindx,paddata(newindx),1:numel(paddata));
        outdata(cTrial,cChan,:) = spline_erp(npadding+1:end-npadding);
    end
end
FT_EEG_SPLINE = FT_EEG;
FT_EEG_SPLINE.dimord = 'rpt_chan_time'; % was 'chan_time_rpt' but should be 'rpt_chan_time' for consistency with conversion to raw format in ft_checkdata
FT_EEG_SPLINE.(field) = outdata;
if isfield(FT_EEG_SPLINE,'sampleinfo');
    FT_EEG_SPLINE = rmfield(FT_EEG_SPLINE,'sampleinfo');
end
if isfield(FT_EEG_SPLINE,'cfg');
    FT_EEG_SPLINE = rmfield(FT_EEG_SPLINE,'cfg');
end