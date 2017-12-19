function avweightstruct = adam_plot_BDM_weights(stats,cfg)
% plots BDM weights for a specific time range or frequency rannge
% outputs average weights for that point, which can subsequently be used
% for statistical testing weight maps against each other using
% compare_weights.m
% returns avweightstruct which also contains the pStruct field (stats)
%
% J.J.Fahrenfort, VU 2015, 2016, 2017
if nargin<2
    disp('cannot plot graph without some settings, need at least 2 arguments:');
    help plot_MVPA;
    return
end

% setting some defaults
plot_order = [];
plotsubjects = false;

% unpack config
v2struct(cfg);

% where does this come from
if isfield(stats(1),'cfg')
    if isfield(stats(1).cfg,'startdir')
        startdir = stats(1).cfg.startdir;
    end
    if isfield(stats(1).cfg,'folder')
        folder = stats(1).cfg.folder;
    end
end

% determine which conditions to plot
if ~isempty(plot_order)
    for cPlot = 1:numel(plot_order)
        statindex = find(strncmpi(plot_order{cPlot},{stats(:).condname},numel(plot_order{cPlot})));
        if isempty(statindex)
            error('cannot find condition name specified in cfg.plot_order');
        end
        newstats(cPlot) = stats(statindex);
    end
    stats = newstats;
end

% make figure?
if ~plotsubjects
    title_text = regexprep(regexprep(regexprep(folder,'\\','/'),regexprep(startdir,'\\','/'),''),'_',' ');
    if numel(stats) == 1
        title_text = title_text(1:find(title_text=='/',1,'last')-1);
    end
    fh = figure('name',title_text);
    % make sure all figures have the same size regardless of nr of subplots
    % and make sure they fit on the screen
    screensize = get(0,'screensize'); % e.g. 1920 * 1080
    UL = (screensize([3 4])-50)./[numSubplots(numel(stats),2) numSubplots(numel(stats),1)];
    if all(UL>[600 400])
        UL=[600 400]; % take this as default
    end
    po=get(fh,'position');
    po(3:4)=UL.*[numSubplots(numel(stats),2) numSubplots(numel(stats),1)];
    po(1:2) = (screensize(3:4)-po(3:4))/2; po(logical([po(1:2)<100 0 0])) = round(po(logical([po(1:2)<100 0 0]))/4); % position in the center, push further to left / bottom if there is little space on horizontal or vertical axis axis
    set(fh,'position',po);
    set(fh,'color','w');
end

% main routine
%if numel(stats)>1
    for cStats=1:numel(stats)
        disp(['plot ' num2str(cStats)]);
        subplot(numSubplots(numel(stats),1),numSubplots(numel(stats),2),cStats);
        avweightstruct(cStats) = subplot_BDM_weights(stats(cStats).weights(cStats),stats(cStats),cfg);
        %[map, ~, cfg] = subplot_MVPA(stats(cStats),cfg,cStats); % all in the first color
        title(regexprep(stats(cStats).condname,'_',' '),'FontSize',10);
    end
% else
%     %map = subplot_MVPA(stats,cfg);
%     avweightstruct = subplot_BDM_weights(weights(cCond),stats(cCond),cfg);
%     if ~plotsubjects
%         title(regexprep(stats.condname,'_',' '),'FontSize',10);
%     end
% end

% get some values and give a title
v2struct(cfg);
dimord = stats(1).settings.dimord;
title_text = ['time ' regexprep(num2str(unique(timelim)),' +',' - ') ' ms.'];
if strcmp(dimord,'time_frequency')
    title_text = [title_text ', frequency: ' regexprep(num2str(unique(freqlim)),' +',' - ') ' Hz'];
end
%if ~isempty(subjlim)
%    title_text = [title_text ', subjects: ' regexprep(num2str(unique(subjlim)),' +',',')];
%end
set(gcf,'name',title_text,'numbertitle','off');


% subfunction that does plotting
function avweightstruct = subplot_BDM_weights(weights,stats,cfg)
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
% get cfg
v2struct(cfg);
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
        connectivity = get_connected_electrodes({stats.settings.chanlocs(:).labels});
        [ clusterPvals, pStruct ] = cluster_based_permutation(subjweights,0,cfg,settings,[],connectivity);
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
    topoplot_jjf(plot_data,convertlocs(chanlocs,'cart2topo'),'maplimits',weightlim,'style','map','electrodes','ptslabels','plotrad',.6,'nosedir','+Y','emarker2',{elecs,'o','k',10,1}); 
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