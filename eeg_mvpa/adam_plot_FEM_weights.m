function avweightstruct = adam_plot_FEM_weights(cfg,stats)
% ADAM_PLOT_FEM_WEIGHTS generates topographical maps of classifier forward encoding weights over a
% specified time interval. Requires as input the output of adam_compute_group_MVPA.
%
% J.J.Fahrenfort, VU 2016, 2019

% main loop
for c=1:numel(stats)
    avweightstruct(c) = subplot_FEM_weights(stats(c),cfg);
end

% subfunction that does plotting
function avweightstruct = subplot_FEM_weights(stats,cfg)
%avweightstruct = [];
%pStruct = [];
% get some settings
weights = stats.weights;
v2struct(weights);
% now unpack settings
nconds = 2;
freqs = 0;
settings = stats.settings;
v2struct(settings);
if iscell(chanlocs)
    chanlocs = chanlocs{1};
end

% setting some graph defaults
subjlim = [];
timelim = 250; % indicate time or period for which to compute weights in ms, if empty plots time with highest average weight
freqlim = []; % indicate freq or freq range for which to compute weights in ms, if empty plots frequency with heighest average weight
weightlim = []; % indicate color range to plot
normalized = true; % can be false (default) or true (z-scoring is performed)
imgtype = [];
indiv_pval = .05;
cluster_pval = .05;
iterations = 1000;
one_two_tailed = 'two';
mpcompcor_method = 'uncorrected';
% then unpack graphsettings too
v2struct(cfg);
pval(1) = indiv_pval;
pval(2) = cluster_pval;

% first a hack to change make sure that whatever is in time is expressed as ms
if mean(times{1}<10)
    times = round(times{1} * 1000);
end

data = indivWeights;

% data: subj (* frequency) * time * electrode * channel_response
% select subjects
if ~isempty(subjlim)
    disp(['only plotting subjects: ' vec2str(subjlim)]);
    data = data(subjlim,:,:,:,:);
end

% normalize or not, space is in the 3rd dimension
if normalized
    data = zscore(data,0,3);
end
avWeights = squeeze(mean(data,1)); % average over subjects

% dimord time_time
% avWeights: time * electrode * channel_response
% dimord frequency_time
% avWeights:  frequency * time * electrode * channel response 

% get out some empties, find max absolute weight values
if isempty(timelim) && isempty(freqlim) && strcmp(dimord,'frequency_time')
    [~,freqlim, timelim] = max2d(mean(mean(abs(avWeights),3),4));
end
if isempty(timelim) && strcmp(dimord,'frequency_time')
    if numel(freqlim) == 1
        freqlim(2) = freqlim(1);
    end
    lowIndex = nearest(freqs,freqlim(1));
    highIndex = nearest(freqs,freqlim(2));
    [~,timelim] = max(mean(mean(mean(abs(avWeights(lowIndex:highIndex,:,:,:)),1),2),4));
    timelim = times(timelim);
elseif isempty(timelim)
    [~,timelim] = max(mean(mean(abs(avWeights),2),3));
    timelim = times(timelim);
end
if isempty(freqlim) && strcmp(dimord,'frequency_time')
    if numel(timelim) == 1
        timelim(2) = timelim(1);
    end
    lowIndex = nearest(times,timelim(1));
    highIndex = nearest(times,timelim(2));
    [~,freqlim] = max(mean(mean(mean(abs(avWeights(:,lowIndex:highIndex,:,:)),2),3),4));
    freqlim = freqs(freqlim);
end
if numel(timelim) == 1
    timelim(2) = timelim(1);
end
if numel(freqlim) == 1
    freqlim(2) = freqlim(1);
end

% extract relevant info
lowIndex = nearest(times,timelim(1));
highIndex = nearest(times,timelim(2));
if strcmp(dimord,'time_time')
    subjweights = squeeze(mean(data(:,lowIndex:highIndex,:,:),2));
    avWeights = squeeze(mean(avWeights(lowIndex:highIndex,:,:),1));
else
    lowfreqIndex = nearest(freqs,freqlim(1));
    highfreqIndex = nearest(freqs,freqlim(2));
    subjweights = squeeze(mean(mean(data(:,lowfreqIndex:highfreqIndex,lowIndex:highIndex,:,:),2),3));
    avWeights = squeeze(mean(mean(avWeights(lowfreqIndex:highfreqIndex,lowIndex:highIndex,:,:),1),2));
end
if isempty(weightlim)
    maxweight = max(abs([min2d(avWeights) max2d(avWeights)]));
    weightlim = [-maxweight maxweight];
elseif numel(weightlim) == 1
    weightlim = [-weightlim weightlim];
end

% figure title
nCond = size(avWeights,2);
title_text = [ regexprep(stats.condname,'_',' ') ': time ' regexprep(num2str(unique(timelim)),' +',' - ') ' ms.'];
if strcmp(dimord,'time_frequency')
    title_text = [title_text ', frequency: ' regexprep(num2str(unique(freqlim)),' +',' - ') ' Hz'];
end
if ~isempty(subjlim)
    title_text = [title_text ', subjects: ' regexprep(num2str(unique(subjlim)),' +',',')];
end
fh = figure;
set(fh, 'Position', get(0,'Screensize'));
set(fh,'color','w');
set(gcf,'name',title_text,'numbertitle','off');

% plot result
for cCond = 1:nCond
    
    % do some statistics
    if strcmpi(mpcompcor_method, 'cluster_based')
        connectivity = get_connected_electrodes({chanlocs(:).labels});
        [clusterPvals, pStruct(cCond)] = cluster_based_permutation(subjweights(:,:,cCond),0,cfg,settings,[],connectivity);
        % pStruct(cCond) = tempStruct;
    elseif strcmpi(mpcompcor_method, 'uncorrected')
        if strcmpi(one_two_tailed,'two')
            [~,clusterPvals] = ttest(subjweights(:,:,cCond),0,'tail','both');
        else
            [~,clusterPvals] = ttest(subjweights(:,:,cCond),0,'tail','right');
        end
    elseif strcmpi(mpcompcor_method, 'fdr')
        disp('fdr not yet implemented here');
        clusterPvals = [];
    else
        clusterPvals = [];
    end
    
    plot_data = avWeights(:,cCond);
    elecs = find(clusterPvals<cluster_pval);
    
    nPlots = nCond;
%     if nPlots <= 8
%         nPlots = 8;
%     end 
    subplot(numSubplots(nPlots,1),numSubplots(nPlots,2),cCond);
    title_text = ['condition ' num2str(cCond)];
    
    % plotting weights
    set(gcf,'color','w');
    title(title_text,'FontSize', 32);
    if strcmp(imgtype,'vec')
        topoplot_jjf(plot_data,chanlocs','maplimits',weightlim,'style','blank','electrodes','on','plotrad',.65,'nosedir','-Y','emarker',{'.','k',20,1},'emarker2',{elecs,'o','k',10,1}); %
    elseif strcmp(imgtype,'png')
        topoplot_jjf(plot_data,chanlocs','maplimits',weightlim,'style','map','electrodes','off','plotrad',.65,'hcolor','none','shading','interp','nosedir','-Y'); %
    else
        topoplot_jjf(plot_data,chanlocs','maplimits',weightlim,'style','map','electrodes','ptslabels','plotrad',.65,'emarker2',{elecs,'o','k',10,1}); % ,'nosedir','-Y'
        cbar('vert');
    end
end

% colormap
cmap  = brewermap([],'*RdBu');
colormap(gcf,cmap);

% return the averages
avweightstruct.chanlocs = chanlocs;
avweightstruct.avWeights = avWeights;
avweightstruct.indivWeights = subjweights;
avweightstruct.timelim = timelim;
avweightstruct.freqlim = freqlim;

if ~exist('avweightstruct','var')
    avweightstruct = [];
end
if ~exist('pStruct','var')
    pStruct = [];
end
avweightstruct.pStruct = pStruct;

if isempty(weightlim)
    sameaxes('xyzc',gcf());
end

