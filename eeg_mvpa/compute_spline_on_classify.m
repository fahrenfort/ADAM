function outdata = compute_spline_on_classify(classify_data,time,splinefreq,relevant_window_in_ms)
% turn classification accuracy into spline
% use 8hz as default frequency
% relevant_window_in_ms sets the window to detect the peak of the spline in
% milliseconds (defaults to [50,750]).
% The spline is fitted around the peak to prevent dilution of peak accuracy
% Johannes Fahrenfort, VU, 2016
samplefreq = numel(time)/max(time) - min(time);
npadding = round(samplefreq/4); % 250 ms
if nargin < 4 || isempty(relevant_window_in_ms)
    relevant_window_in_ms = [50 750];
end

data = classify_data;
outdata = zeros(size(data));
window_indx = nearest(time,relevant_window_in_ms);

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
if ~isempty(splinefreq)
    newstepsize = samplefreq/splinefreq;
else
    newstepsize = 1;
end
% get new indexes based on this sampling frequency
indx = round(1:newstepsize:numel(paddata));
peakdif = indx(nearest(indx,peakindx))-peakindx; % center around max (or min) of trial to get the peak out
newindx = indx - peakdif;
newindx = newindx(newindx>0);
newindx = newindx(newindx<=numel(paddata));
if size(newindx) == 0
    error('this cannot happen');
end
splined_data = spline(newindx,paddata(newindx),1:numel(paddata));
outdata = splined_data(npadding+1:end-npadding);

