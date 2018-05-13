function pstruct = compute_pstruct(labels,clusterPvals,data,cfg,settings,mask,connectivity)
% function pstruct = compute_pstruct(labels,clusterPvals,data,cfg,settings,mask,connectivity)
% computes an array of structs specifying for each cluster:
% - cluster size (number of points in cluster)
% - start_time and stop_time in ms
% - start_freq and stop_freq
% - included electrodes
% - cluster p-value
reduce_dims = '';
v2struct(cfg);
v2struct(settings);
pstruct = [];
if nargin<7
    connectivity = [];
end
if nargin<6
    mask = ones(size(data));
end

mask = logical(mask);
clusterlist = setdiff(unique(labels),0);
for c = 1:numel(clusterlist)
    thisclust = false(size(mask)); % make sure it has the right dimensions
    thisclust(labels == clusterlist(c)) = true;
    pvals = clusterPvals(thisclust);
    pstruct(c).clusterpval = pvals(1);
    pstruct(c).clustersize = sum(sum(thisclust));
    pstruct(c).datasize = sum(sum(mask));
    thisdata = data; thisdata(~thisclust|~mask) = 0; % isolate cluster
    [~,indx1,indx2]=max2d(thisdata); % get peak
    if isempty(connectivity)
        if any(size(data)==1) % this is a 2D plot
            if strcmpi(reduce_dims,'avtrain') || strcmpi(reduce_dims,'avtrain')
                timInd = 2;
            else
                timInd = 1;
            end
            pstruct(c).start_time = round(times{timInd}(find(thisclust(mask), 1,'first'))*1000);
            pstruct(c).stop_time = round(times{timInd}(find(thisclust(mask), 1,'last'))*1000);
            pstruct(c).peak_time = round(times{timInd}(indx1)*1000);
        elseif strcmpi(dimord,'time_time')
            % joram 15-5-17: found a bug here when train and test do not have same
            % time window; the indx1 for one of the clusters falls outside
            % train time window, so it crashes for peak_train; I added a
            % try-catch statement to circumvent the problem
            pstruct(c).start_train = round(times{1}(find(mean(thisclust,1),1,'first'))*1000); % average over test
            pstruct(c).stop_train = round(times{1}(find(mean(thisclust,1),1,'last'))*1000);
            try
                pstruct(c).peak_train = round(times{1}(indx2)*1000);
            catch me
                pstruct(c).peak_train = 0;
            end
            pstruct(c).start_test = round(times{2}(find(mean(thisclust,2),1,'first'))*1000); % average over train
            pstruct(c).stop_test = round(times{2}(find(mean(thisclust,2),1,'last'))*1000);
            try
                pstruct(c).peak_test = round(times{2}(indx1)*1000);
            catch me
                pstruct(c).peak_test = 0;
            end
        elseif strcmpi(dimord,'freq_time')
            freqs = settings.freqs;
            pstruct(c).start_time = round(times{1}(find(mean(thisclust,1),1,'first'))*1000); % average over freq
            pstruct(c).stop_time = round(times{1}(find(mean(thisclust,1),1,'last'))*1000);
            pstruct(c).peak_time = round(times{1}(indx2)*1000);
            if ~strcmpi(reduce_dims,'avfreq')
                pstruct(c).start_freq = freqs(find(mean(thisclust,2),1,'first'));
                pstruct(c).stop_freq = freqs(find(mean(thisclust,2),1,'last'));
                pstruct(c).peak_freq = freqs(indx1);
            end
        end
    else
        % quick hack to take the channel set from the training data
%         if numel(channels) == 2
%             channels = channels{1};
%         end
        % joram 7-2-18: found a bug here, where "channels" has a different label order than
        % chanlocs.labels; the "labels" variable containing the marked channels belonging to a
        % significant cluster are based on the latter, so finding their labels from the "channels"
        % variable resulted in incorrect channel labels in the pStruct (e.g. P7,P5 and PO7 were
        % marked as FC2, FCz and CP2)
        pstruct(c).channels = {chanlocs(logical(labels)).labels};
    end
end