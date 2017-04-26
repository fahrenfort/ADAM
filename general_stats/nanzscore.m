function outdata = nanzscore(data)
% compute zscore on series that may contain NaNs
outdata = (data-nanmean(data))/nanstd(data);