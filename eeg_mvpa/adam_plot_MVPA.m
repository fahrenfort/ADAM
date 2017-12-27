function map = adam_plot_MVPA(cfg,stats)
% function adam_plot_MVPA(cfg,stats)
% draw classification accuracy based on stats produced by
% compute_group_MVPA
%
% By J.J.Fahrenfort, VU, 2014, 2015, 2016, 2017
if nargin<2
    disp('cannot plot graph without stats input, need at least 2 arguments:');
    help plot_MVPA;
    return
end

% settting some defaults
plot_model = 'BDM';
ndec = 2;
inverty = [];
plotsubjects = false;
singleplot = false;
plot_order = [];
folder = '';
startdir = '';
swapaxes = true;
referenceline = 0;   % you can plot a reference line by indicating cfg.referenceline = timepoint (in milliseconds) where it should be plotted.
if numel(stats) > 1
    line_colors = {[.5 0 0], [0 .5 0] [0 0 .5] [.5 .5 0] [0 .5 .5] [.5 0 .5]};
else
    line_colors = {[0 0 0]};
end
acclim3D = [];
acclim2D = [];
cent_acctick = [];

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

% set color-limits (z-axis) or y-limits
if isempty(cent_acctick)
    % assuming this is the same for all graphs, cannot really be helped
    if strcmpi(stats(1).settings.measuremethod,'hr-far') || strcmpi(stats(1).settings.measuremethod,'\muV') || strcmpi(plot_model,'FEM')
        cent_acctick = 0;
    else
        cent_acctick = 1/stats(1).settings.nconds;
    end
end
chance = cent_acctick;

% direction of y-axis
if isempty(inverty)
    if strcmpi(stats(1).settings.measuremethod,'\muV')
        inverty = true;
    else
        inverty = false;
    end
end

% set plottype
if any(size(stats(1).ClassOverTime)==1)
    plottype = '2D'; %  used internally in this function
else
    plottype = '3D';
end
if strcmpi(plottype,'3D')
    singleplot = false;
end

% set line colors
if numel(line_colors)<numel(stats) || isempty(line_colors)
    if numel(stats) > 1
        line_colors = {[.5 0 0], [0 .5 0] [0 0 .5] [.5 .5 0] [0 .5 .5] [.5 0 .5]};
    else
        line_colors = {[0 0 0]};
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

% determine accuracy limits
eval(['acclim = acclim' plottype ';']);
if isempty(acclim)
    mx = max(max(([stats(:).ClassOverTime])));
    mn = min(min(([stats(:).ClassOverTime])));
    if strcmpi(plottype,'2D') % this is a 2D plot
        shift = abs(diff([mn mx]))/10;
        acclim = [mn-shift mx+shift];
    else
        shift = abs(diff([mn mx]))/20;
        mx = max(abs([mx-chance chance-mn]));
        acclim = [-mx-shift mx+shift]+chance;
    end
end

% pack config with defaults
nameOfStruct2Update = 'cfg';
cfg = v2struct(inverty,acclim,chance,cent_acctick,line_colors,ndec,plottype,singleplot,swapaxes,referenceline,nameOfStruct2Update);

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
    if singleplot
        UL = screensize([3 4])-50; % subtract a little to be on the save side
    else
        UL = (screensize([3 4])-50)./[numSubplots(numel(stats),2) numSubplots(numel(stats),1)];
    end
    if all(UL>[600 400])
        UL=[600 400]; % take this as default
    end
    po=get(fh,'position');
    if singleplot
        po(3:4)=UL;
    else
        po(3:4)=UL.*[numSubplots(numel(stats),2) numSubplots(numel(stats),1)];
    end
    po(1:2) = (screensize(3:4)-po(3:4))/2; po(logical([po(1:2)<100 0 0])) = round(po(logical([po(1:2)<100 0 0]))/4); % position in the center, push further to left / bottom if there is little space on horizontal or vertical axis axis
    set(fh,'position',po);
    set(fh,'color','w');
end

% main routine
if numel(stats)>1
    for cStats=1:numel(stats)
        disp(['plot ' num2str(cStats)]);
        if singleplot
            hold on;
            [map, H, cfg] = subplot_MVPA(stats(cStats),cfg,cStats);
            legend_handle(cStats) = H.mainLine;
            legend_text{cStats} = regexprep(stats(cStats).condname,'_',' ');
        else
            subplot(numSubplots(numel(stats),1),numSubplots(numel(stats),2),cStats);
            [map, ~, cfg] = subplot_MVPA(stats(cStats),cfg,cStats); % all in the first color
            title(regexprep(stats(cStats).condname,'_',' '),'FontSize',10);
        end
    end
    if singleplot
        legend(legend_handle,legend_text,'FontSize',10);
        legend boxoff;
    end
else
    map = subplot_MVPA(stats,cfg);
    if ~plotsubjects
        title(regexprep(stats.condname,'_',' '),'FontSize',10);
    end
end

function [map, H, cfg] = subplot_MVPA(stats,cfg,cGraph)
map = [];
if nargin<3
    cGraph = 1;
end

% setting some graph defaults
trainlim = [];
testlim = [];
freqlim = [];
freqtick = [];
plot_model = [];
reduce_dims = [];
inverty = false;
downsamplefactor = 1;
smoothfactor = 1;
plotsigline_method = 'both';
splinefreq = [];
timetick = [];
acctick = [];
mpcompcor_method = 'uncorrected';
plotsubjects = false;
cluster_pval = .05;
indiv_pval = .05;
one_two_tailed = 'two';

% then unpack cfgs to overwrite defaults, first original from stats, then the new one 
if isfield(stats,'cfg')
    oldcfg = stats.cfg;
    v2struct(oldcfg); % unpack the stats-specific cfg
    stats = rmfield(stats,'cfg');
end
v2struct(cfg);
if exist('plot_dim','var') && strcmpi(plot_dim, 'freq_time')
    swapaxes = false;
else
    swapaxes = true;
end
if isempty(splinefreq)
    makespline = false;
else
    makespline = true;
end

% unpack stats
v2struct(stats);
settings = stats.settings; % little hardcoded hack because settings.m is apparently also a function

% now settings should be here, unpack these too
freqs = 0;
v2struct(settings);

% fill some empties
if isempty(freqtick)
    if max(freqlim) >= 60
        freqtick = 20;
    else
        freqtick = 10;
    end
end

% first a hack to change make sure that whatever is in time is expressed as ms
if mean(times{1}<10)
    times{1} = round(times{1} * 1000);
    if numel(times) > 1
        times{2} = round(times{2} * 1000);
    end
end

% determine axes
if strcmpi(reduce_dims,'avtrain') 
    xaxis = times{2};
elseif strcmpi(reduce_dims,'avtest') || strcmpi(reduce_dims,'avfreq')
    xaxis = times{1};
elseif strcmpi(dimord,'time_time')
    xaxis = times{1};
    yaxis = times{2};
elseif strcmpi(dimord,'freq_time')
    xaxis = times{1};
    yaxis=round(freqs);
end

% determine time tick
if isempty(timetick)
    timetickoptions = [5,10,25,50,100,250,500,750,1000,1250,1500,2000,2500,5000];
    timetick = timetickoptions(nearest(timetickoptions,(max(xaxis) - min(xaxis))/5));
end

% if time tick is smaller than sample duration, increase time tick to sample duration
if numel(xaxis)>1 && timetick < (xaxis(2)-xaxis(1))
    timetick = ceil((xaxis(2)-xaxis(1)));
end

% get x-axis and y-axis values
[dims] = regexp(dimord, '_', 'split');
ydim = dims{1};
xdim = dims{2}; % unused

if numel(acclim) == 1
    acclim = [-acclim acclim];
end
acclim = sort(acclim);

% determine acctick
if isempty(acctick)
    if plotsubjects
        mx = max(max((ClassOverTime(:))));
        mn = min(min((ClassOverTime(:))));
        if strcmpi(plottype,'3D')
            acctick = abs(max(abs([chance-mn mx-chance])))-0.001;
        else
            acctick = abs(min(abs([chance-mn mx-chance])));
        end
        if acctick ==0
            acctick = .5;
        end
    else
        if strcmpi(measuremethod,'\muV')
            acctick = round(diff(acclim)/8);
            if acctick ==0
                acctick = 1;
            end
        else
            acctick = diff(acclim)/8;
        end
    end
end

% stuff particular to 2D and 3D plotting
if strcmpi(plottype,'2D')
    yaxis = sort(unique([cent_acctick:-acctick:min(acclim) cent_acctick:acctick:max(acclim)]));
    data = ClassOverTime;
    stdData = StdError;
    % some downsampling on 2D
    if ~isempty(downsamplefactor) && downsamplefactor > 0
        xaxis = downsample(xaxis,downsamplefactor);
        data = downsample(data,downsamplefactor);
        if ~isempty(StdError)
            stdData = downsample(stdData,downsamplefactor);
        end
        if ~isempty(pVals)
            pVals = downsample(pVals,downsamplefactor);
        end
    end
    % some smoothing on 2D using spline
    if makespline
       data = compute_spline_on_classify(data,xaxis',splinefreq);
       if ~isempty(stdData)
           stdData = compute_spline_on_classify(stdData,xaxis',splinefreq);
       end
    end
else
    zaxis = sort(unique([cent_acctick:-acctick:min(acclim) cent_acctick:acctick:max(acclim)]));
    data = ClassOverTime;
end

% set ticks for x-axis and y-axis
xticks = timetick;
if strcmpi(plottype,'3D')
    if strcmpi(ydim,'freq')
        yticks = freqtick;
    else
        yticks = xticks;
    end
else
    yticks = acctick;
end

% some smoothing on 3D, this may neeed fixing?
if makespline && strcmpi(plottype,'3D')
    xaxis = spline(1:numel(xaxis),xaxis,linspace(1, numel(xaxis), round(numel(xaxis)/smoothfactor))); 
    yaxis = spline(1:numel(yaxis),yaxis,linspace(1, numel(yaxis), round(numel(yaxis)/smoothfactor))); 
end

% make a timeline that has 0 as zero-point and makes steps of xticks
if min(xaxis) < 0 && max(xaxis) > 0
    findticks = sort(unique([0:-xticks:min(xaxis) 0:xticks:max(xaxis)]));
else
    findticks = sort(unique(min(xaxis):xticks:max(xaxis)));
end
indx = [];
for tick = findticks
    indx = [indx nearest(xaxis,tick)];
end

% do the same for y-axis
if strcmpi(ydim,'freq')
    findticks = yticks:yticks:max(yaxis);
else
    if min(yaxis) < 0 && max(yaxis) > 0
        findticks = sort(unique([0:-yticks:min(yaxis) 0:yticks:max(yaxis)]));
    else
        findticks = sort(unique(min(yaxis):yticks:max(yaxis)));
    end
end
indy = [];
for tick = findticks
    indy = [indy nearest(yaxis,tick)];
end

% plot
if plotsubjects
    fontsize = 10;
else
    fontsize = 14;
end
if strcmpi(plottype,'2D')
    colormap('default');
    if isempty(StdError)
        H.dataLine = plot(data);
    else
        %H = shadedErrorBar(1:numel(data),data,stdData,{'Color',[.7,.7,.7],'MarkerFaceColor',[1 1 0]},.5); % [1 1 0] {'MarkerFaceColor',[.7,.7,.7]}
        H = shadedErrorBar(1:numel(data),data,stdData,{'Color',line_colors{cGraph} + (1 - line_colors{cGraph})/2,'MarkerFaceColor',[1 1 0]},.1);
        %H.dataLine = plot(data,'Color',line_colors{cGraph}/3,'LineWidth',1);
    end
    hold on;
    % plot horizontal line on zero
    plot([1,numel(data)],[chance,chance],'k--');
    % plot vertical line on zero
    plot([nearest(xaxis,0),nearest(xaxis,0)],[acclim(1),acclim(2)],'k--');
    
    % plot help line
    if referenceline ~= 0
        plot([nearest(xaxis,referenceline),nearest(xaxis,referenceline)],[acclim(1),acclim(2)],'k--');
    end
    
    % plot significant time points
    if ~isempty(pVals)
        sigdata = data;
        if strcmpi(plotsigline_method,'straight') || strcmpi(plotsigline_method,'both') && ~plotsubjects && ~strcmpi(mpcompcor_method,'none')
            if ~singleplot elevate = 1; else elevate = cGraph-.5; end
            if inverty
                sigdata(1:numel(sigdata)) = max(acclim) - (diff(acclim)/80)*elevate;
            else
                sigdata(1:numel(sigdata)) = min(acclim) + (diff(acclim)/80)*elevate;
            end
            sigdata(pVals>=indiv_pval) = NaN;
            if isnumeric(line_colors{cGraph})
                H.bottomLine=plot(1:numel(sigdata),sigdata,'Color',line_colors{cGraph},'LineWidth',3); % sigline below graph
            else
                H.bottomLine=plot(1:numel(sigdata),sigdata,line_colors{cGraph},'LineWidth',2); % sigline below graph
            end
        end
        sigdata = data;
        if ~strcmpi(plotsigline_method,'straight')
            sigdata(pVals>=indiv_pval) = NaN;
            if isnumeric(line_colors{cGraph})
                H.mainLine=plot(1:numel(sigdata),sigdata,'Color',line_colors{cGraph},'LineWidth',3); % sigline on graph
            else
                H.mainLine=plot(1:numel(sigdata),sigdata,line_colors{cGraph},'LineWidth',3); % sigline on graph
            end
        end
        if ~all(isnan((sigdata)))
            wraptext('Due to a bug in the way Matlab exports figures (the ''endcaps'' property in OpenGL is set to''on'' by default), the ''significance lines'' near the time line are not correctly plotted when saving as .eps or .pdf. The workaround is to open these plots in Illustrator, manually select these lines and select ''butt cap'' for these lines (under the ''stroke'' property).');
        end
        if strcmpi(one_two_tailed,'one') one_two_tailed = '1'; else  one_two_tailed = '2'; end
        if ~plotsubjects
            if strcmpi(mpcompcor_method,'uncorrected')
                h_legend = legend(H.mainLine,[' p < ' num2str(indiv_pval) ' (uncorrected, ' one_two_tailed '-sided)']); % ,'Location','SouthEast'
            elseif strcmpi(mpcompcor_method,'cluster_based')
                h_legend = legend(H.mainLine,[' p < ' num2str(cluster_pval) ' (cluster based, ' one_two_tailed '-sided)']);
            elseif strcmpi(mpcompcor_method,'fdr')
                h_legend = legend(H.mainLine,[' p < ' num2str(cluster_pval) ' (FDR, ' one_two_tailed '-sided)']);
            end
            if ~strcmpi(mpcompcor_method,'none')
                legend boxoff;
                set(h_legend,'FontSize',12);
            end
        end
    end
    ylim(acclim);
    xlim([1 numel(data)]);
    % little hack
    if strcmpi(plot_model,'FEM')
        measuremethod = 'CTF slope';
    end
    ylabel(measuremethod,'FontSize',fontsize);
    if strcmpi(reduce_dims,'avtrain')
        xlabel('test time in ms','FontSize',fontsize);
    elseif strcmpi(reduce_dims,'avtest')
        xlabel('train time in ms','FontSize',fontsize);
    else
        xlabel('time in ms','FontSize',fontsize);
    end
    set(gca,'YTick',yaxis);
    Ylabel = regexp(deblank(sprintf(['%0.' num2str(ndec) 'f '],yaxis)),' ','split');
    if chance ~= 0 % create labels containing equal character counts when centered on some non-zero value
        Ylabel((yaxis == chance)) = {'chance'}; % say "chance".
    end
    set(gca,'YTickLabel',Ylabel);
    hold off;
else
    % determine significant time points
    %colormap('jet');
    %cmap  = brewermap([],'*RdBu');
    %colormap(cmap);
    if ~isempty(pVals) && ~strcmpi(mpcompcor_method,'none')
        [data, map] = showstatsTFR(data,pVals,acclim);
    end
    
    % some smoothing on 3D
    if makespline
        if ndims(data) > 2 % DOUBLE CHECK WHETHER THIS IS OK
            [X,Y,Z] = meshgrid(1:size(data,2),1:size(data,1),1:size(data,3));
            [XX,YY,ZZ] = meshgrid(linspace(1,size(data,2),round(size(data,2)/smoothfactor)),linspace(1,size(data,1),round(size(data,1)/smoothfactor)),1:size(data,3));
            data = interp3(X,Y,Z,data,XX,YY,ZZ);
        else
            % this cannot happen because we are in 3D plotting part right?
            [X,Y] = meshgrid(1:size(data,2),1:size(data,1));
            [XX,YY] = meshgrid(linspace(1,size(data,2),round(size(data,2)/smoothfactor)),linspace(1,size(data,1),round(size(data,1)/smoothfactor)));
            data = uint8(interp2(X,Y,double(data),XX,YY));
        end
    end
    if swapaxes
        data = permute(data,[2 1 3]);
    end
    imagesc(data);
    caxis(acclim);
    set(gca,'YDir','normal'); % set the y-axis right
    
    % plot some help lines
    if referenceline ~= 0
        hold on;
        timeinms = referenceline;
        plot([nearest(xaxis,timeinms),nearest(xaxis,timeinms)],[nearest(yaxis,min(yaxis)),nearest(yaxis,max(yaxis))],'k--');
        plot([nearest(xaxis,min(xaxis)),nearest(xaxis,max(xaxis))],[nearest(yaxis,timeinms),nearest(yaxis,timeinms)],'k--');
    end
    
    % set ticks on color bar
    hcb=colorbar;
    Ylabel = regexp(deblank(sprintf(['%0.' num2str(ndec) 'f '],zaxis)),' ','split');
    if chance ~= 0 % create labels containing equal character counts when centered on some non-zero value
        Ylabel((zaxis == chance)) = {'chance'}; % say "chance".
    end
    set(hcb,'YTick',zaxis,'YTickLabel',Ylabel);
    if strcmpi(ydim,'freq')
        ylabel('frequency in Hz','FontSize',fontsize);
        xlabel('time in ms','FontSize',fontsize);
    else
        if swapaxes
            ylabel('training time in ms','FontSize',fontsize);
            xlabel('testing time in ms','FontSize',fontsize);
        else
            ylabel('testing time in ms','FontSize',fontsize);
            xlabel('training time in ms','FontSize',fontsize);
        end
    end
    set(gca,'YTick',indy);
    roundto = yticks;
    set(gca,'YTickLabel',num2cell(int64(round(yaxis(indy)/roundto)*roundto))); % convert to int64 to prevent -0 at 0-axis
end
% set ticks on horizontal axis
set(gca,'XTick',indx);
roundto = xticks;
set(gca,'XTickLabel',num2cell(int64(round(xaxis(indx)/roundto)*roundto)));
set(gca,'FontSize',fontsize);
set(gca,'color','none');
axis square;
if inverty
    set(gca,'YDir','reverse');
end
if plotsubjects && strcmpi(plottype,'2D')
    sameaxes('xyzc',gcf());
end

% invent handle if it does not exist
if ~exist('H')
    H = [];
end
