function avweightstruct = adam_plot_BDM_weights(cfg,stats)
% ADAM_PLOT_BDM_WEIGHTS generates topographical maps of classifier weights or topographical
% class-separability correlation or covariance maps (according to Haufe et al., 2014) over a
% specified time interval and/or frequency range. Requires as input the output of
% adam_compute_group_MVPA.
% 
% By default it performs a cluster-based permutation test over channels/sensors, highlighting
% regions that contributed significantly to the classification performance. Note that this function
% can only be performed if all channels/sensors were selected during adam_MVPA_firstlevel (i.e.
% channelpool = 'ALL' or 'ALL_NOSELECTION'). Besides plotting, this function also returns a stats
% output structure which contains the statistical result of the tested topography. Note that (as the
% name of the function already suggests), this function can only be used for backward decoding and
% not forward encoding models. Also note that the testing time dimension is relevant here (i.e. only
% the diagonal of a time-time generalization analysis), as this concerns the classifier weights and
% not classifier performance.
%
% Use as:
%   adam_plot_BDM_weights(cfg);
%
% The cfg (configuration) input structure can contain the following:
%  
%       cfg.plotweights_or_pattern  = 'weights' (default); this plots the raw classifier weights;
%                                     alternative options are 'covpattern' or 'corpattern'.
%       cfg.timelim                 = 250 (default); or other integer or [int int] range, in ms,
%                                     specifying the desired timepoint or time window to plot
%       cfg.freqlim                 = [int int]; default: [] empty; if you have a time-by-frequency
%                                     diagonal matrix, you can specify a frequency range
%       cfg.mpcompcor_method        = 'cluster_based' (default); string specifying the method for 
%                                     multiple correction correction; other options are:
%                                     'uncorrected', 'fdr' for false-discovery rate, or 'none' if
%                                     you don't wish to perform a statistical analysis. note that
%                                     'cluster_based' makes use of the Fieldtrip toolbox and its
%                                     neighbour/connectivity functions to determine which
%                                     channels/sensors can form clusters.
%       cfg.indiv_pval              = .05 (default); integer; the statistical threshold for each 
%                                     individual channel/sensor; the fdr correction is applied on
%                                     this threshold.
%       cfg.cluster_pval            = .05 (default); integer; if mpcompcor_method is set to
%                                     'cluster_based', this is the statistical threshold for
%                                     evaluating whether a cluster of significant contiguous
%                                     channels/sensors (after the indiv_pval threshold) is larger
%                                     than can be expected by chance; the cluster_pval should never
%                                     be higher than the indiv_pval.
%       cfg.normalized              = true (default); if true, a zscore is performed to have mean 
%                                     across channels/sensors of 0; this facilitates statistical
%                                     analysis of local clusters of classifier weights; can also be
%                                     set to false.
%       cfg.plot_order              = string e.g. {'cond1','cond2'} to specify which conditions you 
%                                     want to extract (and in which order).
%       cfg.nosedir                 = ['+X'|'-X'|'+Y'|'-Y'] direction of nose (default: '+X').
%
% The output weight_stats structure will contain the following fields:
%
%       weight_stats.chanlocs:      [1xN struct]; for N channels/sensors, describing location,
%                                   label, etc.
%       weight_stats.avWeights:     [1xN double]; group-average weights (or covpattern/corpattern)
%       weight_stats.indivWeights:  [MxN double]; subject-specific weigths for M subjects
%       weight_stats.timelim:       [int int]; timelim as specified in input cfg
%       weight_stats.freqlim:       [int int]; freqlim if specified in input cfg
%       weight_stats.pStruct:       [1x1 struct]; p-value statistics
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Example usage: 
%
% cfg = [];
% cfg.plotweights_or_pattern = 'corpattern';
% cfg.timelim                = [100 300];
% cfg.mpcompcor_method       = 'cluster_based';
% cfg.indiv_pval             = .05;
% cfg.cluster_pval           = .01;
%
% weight_stats = adam_plot_BDM_weights(cfg,stats);
%
% part of the ADAM toolbox, by J.J.Fahrenfort, VU, 2017/2018
%
% See also ADAM_COMPUTE_GROUP_MVPA, ADAM_MVPA_FIRSTLEVEL, ADAM_PLOT_MVPA, ADAM_COMPARE_MVPA_STATS

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

% BACKWARDS COMPATIBILITY
if exist('one_two_tailed','var')
    error('The cfg.one_two_tailed field has been replaced by the cfg.tail field. Please replace cfg.one_two_tailed with cfg.tail using ''both'', ''left'' or ''right''. See help for further info.');
end

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
    if all(UL>[400 400])
        UL=[400 400]; % take this as default
    end
    po=get(fh,'position');
    po(3:4)=UL.*[numSubplots(numel(stats),2) numSubplots(numel(stats),1)];
    po(1:2) = (screensize(3:4)-po(3:4))/2; po(logical([po(1:2)<100 0 0])) = round(po(logical([po(1:2)<100 0 0]))/4); % position in the center, push further to left / bottom if there is little space on horizontal or vertical axis axis
    set(fh,'position',po);
    set(fh,'color','w');
end

% main routine
avweightstruct = [];
for cStats=1:numel(stats)
    disp(['plot ' num2str(cStats)]);
    subplot(numSubplots(numel(stats),1),numSubplots(numel(stats),2),cStats);
    avweightstruct = [avweightstruct subplot_BDM_weights(cfg,stats(cStats))];
end


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
function avweightstruct = subplot_BDM_weights(cfg,stats)
weights = stats.weights;
% get some settings
clusterPvals = [];
pStruct = [];
v2struct(weights);
% set defaults
subjlim = [];
timelim = 250;
freqlim = [];
plotweights_or_pattern = 'weights';
normalized = true;
weightlim = 'absmax';
imgtype = [];
indiv_pval = .05;
cluster_pval = .05;
iterations = 1000;
tail = 'both';
mpcompcor_method = 'cluster_based';
nosedir = '+X';
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
        [~,clusterPvals] = ttest(subjweights,0,'tail',tail);
    elseif strcmpi(mpcompcor_method, 'fdr')
        error('fdr not yet implemented here');
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
title(title_text,'FontSize', 12);
if strcmp(imgtype,'vec') % chanlocs' +Y
    topoplot_jjf(plot_data,convertlocs(chanlocs,'cart2topo'),'maplimits',weightlim,'style','blank','electrodes','on','nosedir',nosedir,'emarker',{'.','k',5,1},'emarker2',{elecs,'.','k',15,1}); %
elseif strcmp(imgtype,'png')
    topoplot_jjf(plot_data,convertlocs(chanlocs,'cart2topo'),'maplimits',weightlim,'style','map','electrodes','off','shading','interp','nosedir',nosedir,'hcolor','none');
else
    topoplot_jjf(plot_data,convertlocs(chanlocs,'cart2topo'),'maplimits',weightlim,'style','map','electrodes','on','nosedir',nosedir,'emarker2',{elecs,'o','k',10,1}); % 'electrodes','ptslabels', 'plotrad',.7
    h = cbar('vert');
    set(get(h,'title'),'string',' ');
    %h = colorbar;
    if normalized
        bartitle = 'spatially z-scored';
    else
        bartitle = 'arbitrary units';
    end
    ylabel(h,bartitle);
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
avweightstruct.pStruct = pStruct;
