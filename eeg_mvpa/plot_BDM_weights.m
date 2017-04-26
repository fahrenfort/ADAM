function avweightstruct = plot_BDM_weights(weights,stats,gsettings)
% plots BDM weights for a specific time range or frequency rannge
% outputs average weights for that point, which can subsequently be used
% for statistical testing weight maps against each other using
% compare_weights.m
% returns avweightstruct which also contains the pStruct field (stats)
% J.J.Fahrenfort, VU 2015, 2016

% main loop
fh = figure;
set(fh, 'Position', get(0,'Screensize'));

% plot result
for cCond=1:numel(stats)
    nPlots = numel(stats);
    if nPlots <= 8
        nPlots = 8;
    end
    subplot(numSubplots(nPlots,1),numSubplots(nPlots,2),cCond);
    % plot weights
    avweightstruct(cCond) = subplot_BDM_weights(weights(cCond),stats(cCond),gsettings);
end
% get some values and give a title
subjlim = [];
v2struct(gsettings);
dimord = stats(1).settings.dimord;
title_text = ['time ' regexprep(num2str(unique(timelim)),' +',' - ') ' ms.'];
if strcmp(dimord,'time_frequency')
    title_text = [title_text ', frequency: ' regexprep(num2str(unique(freqlim)),' +',' - ') ' Hz'];
end
if ~isempty(subjlim)
    title_text = [title_text ', subjects: ' regexprep(num2str(unique(subjlim)),' +',',')];
end
set(gcf,'name',title_text,'numbertitle','off');


% subfunction that does plotting
function avweightstruct = subplot_BDM_weights(weights,stats,gsettings)
% get some settings
clusterPvals = [];
pStruct = [];
v2struct(weights);
% set defaults
subjlim = [];
timelim = 250;
freqlim = [];
plotweights_or_pattern = 'corpattern';
normalized = false;
weightlim = 'absmax';
imgtype = [];
indiv_pval = .05;
cluster_pval = .05;
iterations = 1000;
one_two_tailed = 'two';
mpcompcor_method = 'cluster_based';
% get gsettings
v2struct(gsettings);
pval(1) = indiv_pval;
pval(2) = cluster_pval;
% now unpack settings
nconds = 2;
freqs = 0;
settings = stats.settings;
v2struct(settings);

% first a hack to change make sure that whatever is in time is expressed as ms
if mean(times{1}<10)
    times = round(times{1} * 1000);
end

% indivWeights: subj (* frequency) * time * electrode * channel_response
% extract all pattern or weight values
if strcmp(plotweights_or_pattern, 'weight') || strcmp(plotweights_or_pattern, 'weights')
    data = weights.indivWeights;
elseif  strcmp(plotweights_or_pattern, 'covpattern') 
    data = weights.indivCovPatterns;
else
    data = weights.indivCorPatterns;
end

% data: subj (* frequency) * time * electrode * channel_response
% select subjects
if ~isempty(subjlim)
    data = squeeze(data(subjlim,:,:,:,:));
end

% normalize or not, space is in the 3rd dimension
if normalized
    data = zscore(data,0,3);
end
avWeights = squeeze(mean(data,1)); % average over subjects

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

% do some statistics
if isempty(clusterPvals)
    if strcmpi(mpcompcor_method, 'cluster_based')
        connectivity = get_connected_electrodes({weights(1).chanlocs(:).labels});
        [ clusterPvals, pStruct ] = cluster_based_permutation(subjweights,0,gsettings,settings,[],connectivity);
    elseif strcmpi(mpcompcor_method, 'uncorrected')
        if strcmpi(one_two_tailed,'two')
            [~,clusterPvals] = ttest(subjweights,0,'tail','both');
        else
            [~,clusterPvals] = ttest(subjweights,0,'tail','right');
        end
    elseif strcmpi(mpcompcor_method, 'fdr')
        disp('fdr not yet implemented here');
        clusterPvals = [];
    else
        clusterPvals = [];
    end
end

plot_data = avWeights;
elecs = find(clusterPvals<cluster_pval);

% figure title
title_text = {regexprep(stats.condname,'_',' '),['at ' regexprep(num2str(unique(timelim)),' +',' - ') ' ms.']};
if strcmp(dimord,'time_frequency')
    title_text{end+1} = ['frequency: ' regexprep(num2str(unique(freqlim)),' +',' - ') ' Hz'];
end
if ~isempty(subjlim)
    title_text{end+1} = ['subjects: ' regexprep(num2str(unique(subjlim)),' +',',')];
end

% plotting weights
set(gcf,'color','w');
title(title_text,'FontSize', 32);
if strcmp(imgtype,'vec')
    topoplot_jjf(plot_data,chanlocs','maplimits',weightlim,'style','blank','electrodes','on','plotrad',.6,'nosedir','-Y','emarker',{'.','k',20,1},'emarker2',{elecs,'o','k',10,1}); %
elseif strcmp(imgtype,'png')
    topoplot_jjf(plot_data,chanlocs','maplimits',weightlim,'style','map','electrodes','off','plotrad',.6,'hcolor','none','shading','interp','nosedir','-Y'); %
else
    topoplot_jjf(plot_data,chanlocs','maplimits',weightlim,'style','map','electrodes','ptslabels','plotrad',.6,'nosedir','-Y','emarker2',{elecs,'o','k',10,1}); 
    cbar('vert');
end

% return the averages
avweightstruct.chanlocs = chanlocs;
avweightstruct.avWeights = avWeights;
avweightstruct.indivWeights = subjweights;
avweightstruct.timelim = timelim;
avweightstruct.freqlim = freqlim;
avweightstruct.pStruct = pStruct;

%if isempty(weightlim)
%    sameaxes('xyzc',gcf());
%end
