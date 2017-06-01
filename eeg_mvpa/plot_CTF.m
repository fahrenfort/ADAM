function avgctfstruct = plot_CTF(stats,cfg)
% function plot_CTF(stats,weights,cfg)
% plot channel responses of forward encoding model, takes data from
% compute_group_MVPA
% THIS FUNCTION NEEDS SOME UPDATING TO RELY ONLY ON compute_group_MVPA TO
% LIMIT DIMENSIONS AND COMPUTE STATS, AND MORE GENERALLY NEEDS A THOROUGH
% CLEANUP
%
% By J.J.Fahrenfort, VU 2016 
if nargin<2
    disp('cannot plot graph without some settings, need at least 2 arguments:');
    help plot_CTF;
    return;
end
timelim = [];
plotfield = 'CTF';
folder = '';
startdir = '';
singleplot = false;
BLtime = [];
reduce_dims = [];
plottype = '3D';
if isfield(cfg,'BLtime') && cfg.BLtime(1)<= cfg.BLtime(2)
    if numel(stats) > 1
        line_colors = {[0.5 0 0],[1 .5 .5],[0 .5 0],[.5 1 .5],[0 0 .5],[.5 .5 1],[.5 .5 0],[1 1 .5],[0 .5 .5],[.5 1 1],[.5 0 .5],[1 .5 1]};
    else
        line_colors = {[0 0 0],[.9 .9 .9]};
    end
else
    if numel(stats) > 1
        line_colors = {[.5 0 0], [0 .5 0] [0 0 .5] [.5 .5 0] [0 .5 .5] [.5 0 .5]};
    else
        line_colors = {[0 0 0]};
    end
end
v2struct(cfg);
if strcmpi(plottype,'3D')
    singleplot = false;
    cfg.singleplot = false;
    if isempty(reduce_dims) || ~any(strcmpi(reduce_dims,{'diag', 'avtrain', 'avtest'}))
        error('Cannot plot a 3D CTF without averaging over a dimension, specify cfg.reduce_dims = ''diag'', ''avtrain'' or ''avtest'' (use cfg.plottype = ''2D'' to use ''avtraintest'')');
    end
end
if isempty(reduce_dims) || ~any(strcmpi(reduce_dims,{'diag', 'avtrain', 'avtest', 'avtraintest'}))
    error('Cannot plot CTF without averaging over a dimension, specify cfg.reduce_dims = ''diag'', ''avtrain'', ''avtest'' or ''avtraintest''');
end
containsbaseline = ~isempty(BLtime) && cfg.BLtime(1)<= cfg.BLtime(2);
if containsbaseline && numel(line_colors)<numel(stats)*2 || isempty(line_colors)
    if numel(stats) > 1
        line_colors = {[0.5 0 0],[1 .5 .5],[0 .5 0],[.5 1 .5],[0 0 .5],[.5 .5 1],[.5 .5 0],[1 1 .5],[0 .5 .5],[.5 1 1],[.5 0 .5],[1 .5 1]};
    else
        line_colors = {[0 0 0],[.5 .5 .5]};
    end
elseif ~containsbaseline && numel(line_colors)<numel(stats) || isempty(line_colors)
    if numel(stats) > 1
        line_colors = {[.5 0 0], [0 .5 0] [0 0 .5] [.5 .5 0] [0 .5 .5] [.5 0 .5]};
    else
        line_colors = {[0 0 0]};
    end
end
cfg.line_colors = line_colors;

v2struct(stats(1).settings,{'fieldNames','dimord'});

    
% can only use diag if testlim is trainlim and dimension is time_time
if strcmpi(reduce_dims, 'diag') && strcmpi(dimord,'time_time')
    if numel(testlim) == 1
        cfg.testlim = [];
    end
    if numel(trainlim) == 1
        cfg.trainlim = [];
    end
    cfg.testlim = cfg.trainlim;
end
% general time limit
if ~isempty(timelim) % timelim takes precedence
    cfg.trainlim = timelim;
    cfg.testlim = timelim;
end
% limit freq in this case CTF = freq * time * channelweight becomes CTF = time * channelweight
if strcmp(dimord,'freq_time')
    if isempty(freqlim)
        freqlim = input('What frequency or frequency range should I extract? ');
    end
    if numel(freqlim) == 1
        freqlim(2) = freqlim(1);
    end
    cfg.freqlim = freqlim;
end

if strcmp(plotfield,'CTF')
    % make figure
    title_text = regexprep(regexprep(folder,startdir,''),'_',' ');
    title_text = title_text(1:find(title_text==filesep,1,'last')-1);
    fh = figure('name',title_text);
    UL=[600 450];
    po=get(fh,'position');
    % the line below needs to be adjusted for singleplot
    if singleplot
        po(3:4)=UL;
    else
        po(3:4)=UL.*[numSubplots(numel(stats),2) numSubplots(numel(stats),1)];
    end
    set(fh,'position',po);
    set(fh,'color','w');
end

% if singleplot
%     colororder = get(gca,'ColorOrder');
% end

% loop for main conditions
for cStats=1:numel(stats)
    if strcmp(cfg.plotfield,'CTF')
        if singleplot
            hold on;
            if ~isempty(BLtime) && BLtime(1)<= BLtime(2)
                legend_text{cStats*2-1} = ['CTF ' strrep(stats(cStats).condname,'_',' ')];
                legend_text{cStats*2} = 'baseline';
            else
                legend_text{cStats} = ['CTF ' strrep(stats(cStats).condname,'_',' ')];
            end
        else
            if ~isempty(BLtime) && BLtime(1)<= BLtime(2)
                legend_text{1} = 'CTF';
                legend_text{2} = 'baseline';
            else
                legend_text{1} = 'CTF';
            end
            subplot(numSubplots(numel(stats),1),numSubplots(numel(stats),2),cStats);
        end
    end
    [indivCTFs, indivCTFmean, indivBLmean, pStruct] = plot_routine(stats(cStats),cfg,cStats);
    avgctfstruct(cStats).indivCTFs = indivCTFs;
    avgctfstruct(cStats).indivCTFmean = indivCTFmean;
     avgctfstruct(cStats).indivBLmean = indivBLmean;
    avgctfstruct(cStats).pStruct = pStruct;
    avgctfstruct(cStats).condname = stats(cStats).condname;
end

if singleplot
    legend(legend_text);
    legend boxoff;
end

% what to plot: individual condition CTFs or average CTF?
function [indivCTFs, indivCTFmean, indivBLmean, pStruct] = plot_routine(stats,cfg,cGraph)
if nargin<3
    cGraph = 1;
end
indivCTFs = [];
plotfield = 'CTF';
shiftindiv = false;
v2struct(cfg);
weights = stats.weights;
nCond = size(weights.CTF,ndims(weights.CTF));
if strcmp(plotfield,'CTFpercond')
    v2struct(weights);
    % make figure
    fh = figure;
    set(fh, 'Position', get(0,'Screensize'));
    set(fh,'color','w');
    % loop for individual conditions (channels)
    for cCond=1:nCond
        % shift shape of CTF to central point, a bit complicated due to
        % the fact that the array is time*time*chan (or freq*time*chan) which has to be
        % shifted to be chan*time*time (or chan*freq*time) prior to shifting
        channel_pos{cCond} =  1:nCond;
        if shiftindiv
            dim2shift = ndims(CTFpercond{cCond}); % channels are always the last dimension
            CTFpercond{cCond} = circshift(CTFpercond{cCond},floor(nCond)/2-cCond,dim2shift); 
            semCTFpercond{cCond} = circshift(semCTFpercond{cCond},floor(nCond)/2-cCond,dim2shift);
            channel_pos{cCond} = circshift(channel_pos{cCond},floor(nCond)/2-cCond,2);
        end
        weights.CTF = CTFpercond{cCond};
        weights.semCTF = semCTFpercond{cCond};
        weights.channel_pos = channel_pos{cCond};
        subplot(numSubplots(nCond,1),numSubplots(nCond,2),cCond);
        [indivCTFs{cCond}, indivCTFmean{cCond}, indivBLmean{cCond}, pStruct{cCond}] = subplot_CTF(stats,weights,cfg);
        set(gcf,'name',stats.condname,'numbertitle','off');
        title_text = ['condition ' num2str(cCond)];
        title(title_text);
    end
else
    weights.channel_pos =  1:nCond;
    [indivCTFs, indivCTFmean, indivBLmean, pStruct] = subplot_CTF(stats,weights,cfg,cGraph);
    if ~singleplot % ~isfield(cfg,'color')
        title(strrep(stats.condname,'_',' '));
    end
end


% use subfunction to do all the plotting, plots individual CTFs
function [indivCTF, indivCTFmean, indivBLmean, pStruct] = subplot_CTF(stats,weights,cfg,cGraph)
if nargin<4
    cGraph = 1;
end
indivCTFmean = [];
indivBLmean = [];
pStruct = [];
% unpack weights
v2struct(weights);
% now unpack settings
nconds = 2;
freqs = 0;
settings = stats.settings;
v2struct(settings);

% setting some graph defaults
temporalsmoothness = [];
trainlim = []; % select time point for train data, if empty take all
testlim = []; % select time point for test data, if empty take all, if both empty take diag
freqlim = []; % if data is frequency data, always select along freqlim
CTFtime = []; % average over this interval for CTF
BLtime = []; % average over this interval for non-CTF baseline
plottype = '2D'; % '2D' (average over time) or '3D' (time line)
timetick = 250;
weightlim = [];
colorlim = [];
reduce_dims = 'diag'; % 'diag', 'avtrain', 'avtest', 'avfreq'
timetick = 250;
makeround = [];
plotsubjects = false;
flat = true;
indiv_pval = .05;
cluster_pval = .05;
iterations = 1000;
one_two_tailed = 'two';
mpcompcor_method = 'cluster_based';
% then unpack graphsettings too
v2struct(cfg);
pval(1) = indiv_pval;
pval(2) = cluster_pval;
if isempty(colorlim)
    colorlim = weightlim;
end
if numel(BLtime) == 1
    BLtime = [BLtime BLtime];
end

% if time tick is smaller than sample duration, increase time tick to sample duration
if numel(times{1})>1 && timetick < (times{1}(2)-times{1}(1))*1000
    timetick = ceil((times{1}(2)-times{1}(1))*1000);
end

% first a hack to change make sure that whatever is in time is expressed as ms
if mean(times{1}<10)
    times{1} = round(times{1} * 1000);
    if numel(times) > 1
        times{2} = round(times{2} * 1000);
    end
end
% enter back into settings
for c = 1:numel(times)
    settings.times{c} = times{c}/1000; % in seconds again
end

% recompute standard error, to be safe
semCTF = squeeze(std(indivCTF))/sqrt(size(indivCTF,1));

% limit subjects, frequency and time
if strcmp(dimord,'freq_time')
    lowIndex = nearest(freqs,freqlim(1));
    highIndex = nearest(freqs,freqlim(2));
    freqs = freqs(lowIndex:highIndex);
    if ~ismatrix(CTF)
        CTF = squeeze(mean(CTF(lowIndex:highIndex,:,:,:),1));
        semCTF = squeeze(mean(semCTF(lowIndex:highIndex,:,:,:),1));
        indivCTF = squeeze(mean(indivCTF(:,lowIndex:highIndex,:,:,:),2));
    end
%     if numel(freqs)>1
%         % limit time
%         if ~isempty(timelim)
%             if numel(timelim) == 1
%                 timelim(2) = timelim(1);
%             end
%             lowIndex = nearest(times{1},timelim(1));
%             highIndex = nearest(times{1},timelim(2));
%             times{1} = times{1}(lowIndex:highIndex);
%             CTF = CTF(:,lowIndex:highIndex,:,:);
%             semCTF = semCTF(:,lowIndex:highIndex,:,:);
%             indivCTF = indivCTF(:,:,lowIndex:highIndex,:,:);
%         end
%     else
        if ~isempty(timelim)
            if numel(timelim) == 1
                timelim(2) = timelim(1);
            end
            lowIndex = nearest(times{1},timelim(1));
            highIndex = nearest(times{1},timelim(2));
            times{1} = times{1}(lowIndex:highIndex);
            settings.times{1} = times{1}/1000; % hack back into settings
            CTF = CTF(lowIndex:highIndex,:,:,:);
            semCTF = semCTF(lowIndex:highIndex,:,:,:);
            indivCTF = indivCTF(:,lowIndex:highIndex,:,:,:);
        end
%     end
else
    % limit y-axis time
    if ~isempty(testlim)
        if numel(testlim) == 1
            testlim(2) = testlim(1);
        end
        lowIndex = nearest(times{2},testlim(1));
        highIndex = nearest(times{2},testlim(2));
        times{2} = times{2}(lowIndex:highIndex);
        settings.times{2} = times{2}/1000; % hack back into settings
        CTF = CTF(:,lowIndex:highIndex,:,:);
        semCTF = semCTF(:,lowIndex:highIndex,:,:);
        indivCTF = indivCTF(:,:,lowIndex:highIndex,:,:);
    end
    % limit x-axis time
    if ~isempty(trainlim)
        if numel(trainlim) == 1
            trainlim(2) = trainlim(1);
        end
        lowIndex = nearest(times{1},trainlim(1));
        highIndex = nearest(times{1},trainlim(2));
        times{1} = times{1}(lowIndex:highIndex);
        settings.times{1} = times{1}/1000; % hack back into settings
        CTF = CTF(lowIndex:highIndex,:,:,:);
        semCTF = semCTF(lowIndex:highIndex,:,:,:);
        indivCTF = indivCTF(:,lowIndex:highIndex,:,:,:);
    end
end


% actual extraction, either diagonal or average over one of the two time lines
xaxis=times{1}; % default
% extract specific train time, test time or frequency
if strcmpi(dimord,'freq_time')
    if numel(freqs) > 1
        CTF = squeeze(mean(CTF,1));
        semCTF = squeeze(mean(semCTF,1));
        indivCTF = squeeze(mean(indivCTF,2));
    end
else
    if strcmpi(reduce_dims, 'diag')
        for c1 =1:size(CTF,3)
            diagCTF(:,c1) = diag(CTF(:,:,c1));
            for csubj = 1:size(indivCTF,1)
                diagindivCTF(csubj,:,c1) = diag(squeeze(indivCTF(csubj,:,:,c1)));
            end
        end
        CTF = diagCTF;
        indivCTF = diagindivCTF;
    elseif strcmpi(reduce_dims, 'avtrain')
        CTF = squeeze(mean(CTF,1));
        indivCTF = squeeze(mean(indivCTF,2));
        xaxis=times{2};
    elseif strcmpi(reduce_dims, 'avtest')
        CTF = squeeze(mean(CTF,2));
        indivCTF = squeeze(mean(indivCTF,3));
    end
end

% make circular CTF in the case of uneven number of channels
if isempty(makeround)
    if mod(size(CTF,2),2) == 0
        makeround = true;
    else
        makeround = false;
    end
end

if makeround
    if ndims(indivCTF) == 3
        indivCTF = cat(3,indivCTF(:,:,end),indivCTF);
        CTF = cat(2,CTF(:,end),CTF);
    elseif ndims(indivCTF) == 4
        indivCTF = cat(4,indivCTF(:,:,:,end),indivCTF);
        CTF = cat(3,CTF(:,:,end),CTF);
    end
    channel_pos = [ channel_pos(end) channel_pos ];
end

% if 2D
if strcmpi(plottype,'2D')
    % fill some empties
    if isempty(CTFtime) || isempty(BLtime)
        if strcmpi(reduce_dims, 'avtrain') || strcmpi(reduce_dims, 'diag') || strcmpi(reduce_dims, 'avfreq') || strcmpi(reduce_dims, 'avtraintest')
            CTFtime = times{1} > 0;
            BLtime = times{1} <= 0;
        else
            CTFtime = times{2} > 0;
            BLtime = times{2} <= 0;            
        end
    else
        if strcmpi(reduce_dims, 'avtrain') || strcmpi(reduce_dims, 'diag') || strcmpi(reduce_dims, 'avfreq') || strcmpi(reduce_dims, 'avtraintest')
            CTFtime = times{1} >= CTFtime(1) & times{1} <= CTFtime(2);
            BLtime = times{1} >= BLtime(1) & times{1} <= BLtime(2);
        else
            CTFtime = times{2} >= CTFtime(1) & times{2} <= CTFtime(2);
            BLtime = times{2} >= BLtime(1) & times{2} <= BLtime(2);
        end
    end
    
    % save CTF as output
    if strcmpi(reduce_dims, 'avtraintest') % THIS DOES NOT WORK YET -> FIX IT!
        indivCTFmean = squeeze(mean(squeeze(mean(indivCTF(:,CTFtime,CTFtime,:),2)),2));
        indivBLmean = squeeze(mean(squeeze(mean(indivCTF(:,BLtime,BLtime,:),2)),2));
        % figure; plot(indivBLmean'); legend(strsplit(num2str(1:20),' '));figure;
        semCTF = std(squeeze(mean(squeeze(mean(indivCTF(:,CTFtime,CTFtime,:),2)),2)))/sqrt(size(indivCTF,1));
        semBL = std(squeeze(mean(squeeze(mean(indivCTF(:,BLtime,BLtime,:),2)),2)))/sqrt(size(indivCTF,1));
    else
        indivCTFmean = squeeze(mean(indivCTF(:,CTFtime,:,:),2)); % original
        indivBLmean = squeeze(mean(indivCTF(:,BLtime,:,:),2)); % original
        semCTF = std(squeeze(mean(indivCTF(:,CTFtime,:),2)))/sqrt(size(indivCTF,1));
        semBL = std(squeeze(mean(indivCTF(:,BLtime,:),2)))/sqrt(size(indivCTF,1));
    end
    
    % what to plot
    if any(BLtime)
        CTF = [ mean(indivCTFmean,1); mean(indivBLmean,1)]';
        semCTF = [ semCTF; semBL]';
    else
        CTF = mean(indivCTFmean,1)';
        semCTF = mean(semCTF,1)';
    end
    
    % plot
    hold on;
    disp(['plot ' num2str(cGraph)]);
    for cL = 1:size(CTF,2)
        if isnumeric(line_colors{2*cGraph-2+cL})
            errorbar(CTF(:,cL),semCTF(:,cL)/2,'Color',line_colors{2*cGraph-2+cL});
        else
            errorbar(CTF(:,cL),semCTF(:,cL)/2,line_colors{2*cGraph-2+cL});
        end
    end
    if ~isempty(weightlim)
        ylim(weightlim);
        % y-axis
        yticks = .05; % hardcoded for now, fix later
        if min(weightlim) < 0 && max(weightlim) > 0
            yaxis = sort(unique([0:-yticks:min(weightlim) 0:yticks:max(weightlim)]));
        else
            yaxis = sort(unique(min(weightlim):yticks:max(weightlim)));
        end
        set(gca,'Ytick',yaxis);
    end
    ylabel('channel response');
    % set x-axis
    indx = round(linspace(1,size(CTF,1),size(CTF,1)));
    set(gca,'Xtick',indx);
    %set(gca,'XTickLabel',num2cell(mod2base([1:size(CTF,1)]-makeround,size(CTF,1)-makeround)));
    set(gca,'XTickLabel',num2cell(channel_pos));
    xlabel('channel');
    % legend
    if any(BLtime) && ~singleplot %~isfield(cfg,'color');
        legend({'CTF', 'baseline'});
    elseif ~singleplot % ~isfield(cfg,'color');
        legend({'CTF'});
    end
    legend boxoff;
else
    % first do some statistics
    sigline = [];
    if isempty(weightlim)
        templim = minmax(reshape(CTF,1,numel(CTF)));
    else
        templim = weightlim;
    end
    if flat
        maxsig = max(templim);
        addy = 4;
    else
        maxsig = min(templim);
        addy = 10;
    end
    maxChan = round(size(indivCTF,3)/2);
    minChan = size(indivCTF,3);
    if strcmpi(mpcompcor_method, 'cluster_based')
        cfg.plottype = '2D'; % little hack to get pStruct right, because we have removed a dimension
        % [ clusterPvals, pStruct ] = cluster_based_permutation(indivCTF(:,:,maxChan),indivCTF(:,:,minChan),cfg,settings);
        % indivCTF(subj,time,chan)
        % compute slopes
        for cSubj = 1:size(indivCTF,1)
            for cTime = 1:size(indivCTF,2)
                [~,~,slope(cSubj,cTime)] = fit_slope(squeeze(indivCTF(cSubj,cTime,:))');
            end
        end
        [ clusterPvals, pStruct ] = cluster_based_permutation(slope,0,cfg,settings);
        cfg.plottype = '3D'; % little hack to get pStruct right
        sigline = nan(size(clusterPvals));
        sigline(clusterPvals < pval(1)) = maxsig;
    elseif strcmpi(mpcompcor_method, 'uncorrected')
        if strcmpi(one_two_tailed,'two')
            [~,clusterPvals] = ttest(indivCTF(:,:,maxChan),indivCTF(:,:,minChan),'tail','both');
        else
            [~,clusterPvals] = ttest(indivCTF(:,:,maxChan),indivCTF(:,:,minChan),'tail','right');
        end
        sigline = nan(size(clusterPvals));
        sigline(clusterPvals < pval(1)) = maxsig;
    elseif strcmpi(mpcompcor_method, 'fdr')
        if strcmpi(one_two_tailed,'two')
            [~,clusterPvals] = ttest(indivCTF(:,:,maxChan),indivCTF(:,:,minChan),'tail','both');
        else
            [~,clusterPvals] = ttest(indivCTF(:,:,maxChan),indivCTF(:,:,minChan),'tail','right');
        end
        thresh = fdr(squeeze(clusterPvals),pval(2));
        sigline = nan(size(clusterPvals));
        sigline(clusterPvals < thresh) = maxsig;
    end
    
    % interpolate y to make smooth: channels
    y = linspace(1,size(CTF,2),size(CTF,2));
    yq = linspace(1,size(CTF,2),size(CTF,2)*10);
    CTFinterp = spline(y,CTF,yq);
    indy = round(linspace(1,size(CTF,2)*10,size(CTF,2)));
    % indy = indy(1+makeround:end);   
    % interpolate x to make smooth: time  
    if ~isempty(temporalsmoothness) && ~strcmp(dimord,'freq_time')
        sf = round(1000*size(CTF,1)/(xaxis(end)-xaxis(1)));
        downfactor = temporalsmoothness/sf; % resample to temporalsmoothness
        x = linspace(1,size(CTF,1),size(CTF,1));
        xq = linspace(1,size(CTF,1),size(CTF,1)*downfactor);
        CTFinterp = spline(x,CTFinterp',xq); % note: the matrix has been transposed
        upfactor = 250/temporalsmoothness; % upsample to 250 Hz
        x = linspace(1,size(CTFinterp,2),size(CTFinterp,2));
        xq = linspace(1,size(CTFinterp,2),size(CTFinterp,2)*upfactor);
        CTFinterp = spline(x,CTFinterp,xq);
        xaxis = linspace(xaxis(1),xaxis(end),numel(xaxis)*(250/sf));
        % also resample significance line
        if ~isempty(sigline)
            x = linspace(1,numel(sigline),numel(sigline));
            xq = linspace(1,numel(sigline),numel(sigline)*(250/sf));
            sigline = interp1(x,sigline,xq,'nearest');
        end
    else
        CTFinterp = CTFinterp';
    end    
    % plot
    if flat
        imagesc(CTFinterp);
        if ~isempty(sigline)
            hold on;
            sigline(~isnan(sigline)) = 2;
            H.mainLine=line(1:numel(sigline),sigline,'LineWidth',8,'Color','black');
            set(gca,'YDir','normal');
        end
        colorbar;
    else
        surf(CTFinterp,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
        if ~isempty(sigline)
            hold on;
            H.mainLine=plot3(1:numel(sigline),ones(size(sigline))*addy,sigline,'k','LineWidth',8);
        end
    end
    if ~all(isnan((sigline)))
        wraptext('Due to a bug in the way Matlab exports figures (the ''endcaps'' property in OpenGL is set to''on'' by default), the ''significance lines'' near the time line are not correctly plotted when saving as .eps or .pdf. The workaround is to open these plots in Illustrator, manually select these lines and select ''butt cap'' for these lines (under the ''stroke'' property).');
    end
    
    % make a timeline that has 0 as zero-point and makes steps of xticks
    xticks = timetick;
    if min(xaxis) <= 0 && max(xaxis) >= 0
        findticks = sort(unique([0:-xticks:min(xaxis) 0:xticks:max(xaxis)]));
    else
        findticks = sort(unique(min(xaxis):xticks:max(xaxis)));
    end
    indx = [];
    for tick = findticks
        indx = [indx nearest(xaxis,tick)];
    end
    roundto = xticks;
    % draw time
    xlabel('time in ms.');
    set(gca,'XTick',indx);
    xlim([indx(1) indx(end)]);
    set(gca,'XTickLabel',num2cell(round(xaxis(indx)/roundto)*roundto));
    % draw channel numbers
    ylabel('channel');
    set(gca,'Ytick',indy);
    set(gca,'YTickLabel',num2cell(channel_pos));
    ylim([1 indy(end)]);
    % z-axis
    ylabel('channel response');
    if ~isempty(weightlim)
        zlim(weightlim);
    end
    if ~isempty(colorlim)
        caxis(colorlim);
    end
    % flat view?
%     if flat
%         view(2);
%         colorbar;
%     end

    % color scheme
    colormap jet;
    % colormap(brewermap([],'RdBu')); 
end
axis square;
set(gca,'FontSize',16);

if isempty(weightlim)
    sameaxes('xyzc',gcf());
end

indivCTF = squeeze(indivCTF);