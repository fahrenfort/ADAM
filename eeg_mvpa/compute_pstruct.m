function pstruct = compute_pstruct(labels,clusterPvals,data,gsettings,settings,mask,connectivity)
% function pstruct = compute_pstruct(labels,clusterPvals,data,gsettings,settings,mask,connectivity)
% computes an array of structs specifying for each cluster:
% cluster size (number of points in cluster)
% start_time and stop_time in ms
% start_freq and stop_freq
% included electrodes
% cluster p-value
v2struct(gsettings);
v2struct(settings);
pstruct = [];
mask = logical(mask);
clusterlist = setdiff(unique(labels),0);
for c = 1:numel(clusterlist)
    thisclust = labels == clusterlist(c);
    pvals = clusterPvals(thisclust);
    pstruct(c).clusterpval = pvals(1);
    pstruct(c).clustersize = sum(sum(thisclust));
    pstruct(c).datasize = sum(sum(mask));
    thisdata = data; thisdata(~thisclust|~mask) = 0; % isolate cluster
    [~,indx1,indx2]=max2d(thisdata); % get peak
    if isempty(connectivity)
        if any(size(data)==1) % this is a 2D plot
            pstruct(c).start_time = round(times{1}(find(thisclust(mask), 1,'first'))*1000);
            pstruct(c).stop_time = round(times{1}(find(thisclust(mask), 1,'last'))*1000);
            pstruct(c).peak_time = round(times{1}(indx1)*1000);
        elseif strcmpi(dimord,'time_time')
            pstruct(c).start_train = round(times{1}(find(mean(thisclust,1),1,'first'))*1000); % average over test
            pstruct(c).stop_train = round(times{1}(find(mean(thisclust,1),1,'last'))*1000);
            pstruct(c).peak_train = round(times{1}(indx1)*1000);
            pstruct(c).start_test = round(times{2}(find(mean(thisclust,2),1,'first'))*1000); % average over train
            pstruct(c).stop_test = round(times{2}(find(mean(thisclust,2),1,'last'))*1000);
            pstruct(c).peak_test = round(times{1}(indx2)*1000);
        elseif strcmpi(dimord,'freq_time')
            freqs = settings.freqs;
            pstruct(c).start_time = round(times{1}(find(mean(thisclust,1),1,'first'))*1000); % average over freq
            pstruct(c).stop_time = round(times{1}(find(mean(thisclust,1),1,'last'))*1000);
            pstruct(c).peak_time = round(times{1}(indx2)*1000);
            if ~strcmpi(reduce_dims,'avfreq')
                pstruct(c).start_freq = freqs(find(mean(thisclust,2),1,'first'));
                pstruct(c).stop_freq = freqs(find(mean(thisclust,2),1,'last'));
                pstruct(c).peak_train = round(times{1}(indx1)*1000);
            end
        end
    else
        pstruct(c).channels = channels(logical(labels))';
    end
end