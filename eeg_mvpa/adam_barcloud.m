function pVals = adam_barcloud(cfg,inputdata)
% function pVals = adam_barcloud(cfg,inputdata) Makes a grouped bar and/or
% raincloud plot with individual inputdata points connected to each other,
% including p-values above each bar combination. Can be useful when
% plotting data from a 2-factor factorial design. inputdata must have the
% following format:
% inputdata(participants, rep_meas1, rep_meas2);
% rep_meas2 must have at most 2 levels if one wants to perform a t-test
% between those levels
% Various optional arguments can be given to adjust colors, axis names, xtick labels
% using the form:
% defcfg.linecolor = [.6 .6 .6];    % to indicate the color of the lines of individual subjects
% defcfg.conditions = {'1', '2'};   % indicate condition names that the bars refer to
% cfg.facecolor{1} = [0 0.46 0.73]; % colors of the bars
% cfg.facecolor{2} = [0.73 0.46 0];
% defcfg.legendlocation = 'NorthEast'; % where to put the legend (see help legend) 
% defcfg.ylim = [];                 % limits on y-axis
% defcfg.xticklabels = [];          % labels on levels for conditions
% defcfg.ytickspacing = .1;         % spacing between ticks on the y-axis
% defcfg.xlabel = '';               % x-axis label
% defcfg.ylabel = '';               % y-axis label
% defcfg.fontsize = 12;             % font size
% defcfg.title = '';                % title of plot
% defcfg.plotP = true;              % perform t-test and plot p-values, p-values are also returned by the function 
% defcfg.tail = 'both';             % indicate that t-test should be double sided, alternatively: 'right' or 'left' (see help ttest) 
% defcfg.square = true;             % plot the axes as square or not
% defcfg.facealpha = .55;           % transparency
% defcfg.graphtype = 'barcloud';    % can also be 'bar_only', 'cloud_only', or 'dif_cloud'
% defcfg.markersize_mean = 10;
% defcfg.markersize_indiv = 6;

% J.J.Fahrenfort, VU 2022, 2023

% defaults
defcfg.conditions = [];
defcfg.legendlocation = 'NorthWest';
defcfg.facecolor{1} = [0 118 186]/255; % colors of the bars
defcfg.facecolor{2} = [181 23 0]/255;
defcfg.xlim = [];
defcfg.ylim = [];
defcfg.xticklabels = [];
defcfg.ytickspacing = .1;
defcfg.xlabel = '';
defcfg.ylabel = '';
defcfg.title = '';
defcfg.fontsize = 12;
defcfg.plotP = true;    % perform t-test and plot p-values
defcfg.tail = 'both';   % whether t-test is two-sided or one-sided (see help ttest) 
defcfg.square = false;   % plot the axes as square
defcfg.plotindiv = true;
defcfg.facealpha = .55;
defcfg.density_granularity = 200;
defcfg.graphtype = 'barcloud';
defcfg.markersize_mean = 10;
if strcmpi(cfg.graphtype, 'bar_only')
    defcfg.markersize_indiv = 5;
    defcfg.linecolor = [.4 .4 .4];
else
    defcfg.markersize_indiv = 3;
    defcfg.linecolor = [.6 .6 .6];
end

% enter default values
vars = fieldnames(defcfg);
for c=1:numel(vars)
    if ~isfield(cfg,vars{c})
        cfg.(vars{c}) = defcfg.(vars{c});
    end
end
if isempty(cfg.xticklabels)
    for c = 1:size(inputdata,2)
        cfg.xticklabels{c} = sprintf('Level%d',c);
    end
end

% make a version with raincloud plots
for c=1:numel(cfg.facecolor)
    colours(c,:) = cfg.facecolor{c};
end
for c1 = 1:size(inputdata,2)
    for c2 = 1:size(inputdata,3)
        data{c1,c2} = squeeze(inputdata(:,c1,c2));
    end
end

% calculate properties of density plots
[n_plots_per_series, n_series] = size(data);
density_granularity = cfg.density_granularity; 
n_bins = repmat(density_granularity, n_plots_per_series, n_series);
bandwidth = [];

% calculate kernel densities
for i = 1:n_plots_per_series
    for j = 1:n_series
        % compute density using 'ksdensity'
        [ks{i, j}, x{i, j}] = ksdensity(data{i, j}, 'NumPoints', n_bins(i, j), 'bandwidth', bandwidth);
        % Define the faces to connect each adjacent f(x) and the corresponding points at y = 0.
        q{i, j}     = (1:n_bins(i, j) - 1)';
        faces{i, j} = [q{i, j}, q{i, j} + 1, q{i, j} + n_bins(i, j) + 1, q{i, j} + n_bins(i, j)]; 
    end
end

% some admin scaling, offsets etc
if strcmpi(cfg.graphtype, 'cloud_only')
    normalize = num2cell(repmat(3.5 * max(max(cellfun(@max, ks))),size(ks)));
else
    normalize = num2cell(repmat(2.2 * max(max(cellfun(@max, ks))),size(ks)));
end
ks = cellfun(@mrdivide,ks,normalize,'UniformOutput',false);
spacing = 1;
if ~mod(j,2) && any(strcmpi(cfg.graphtype,{'cloud_only'}))
    ks_offsets  = [.8:n_plots_per_series] .* spacing;
else
    ks_offsets  = [.7:n_plots_per_series] .* spacing;
end

% calculate patch vertices from kernel density
for i = 1:n_plots_per_series
    for j = 1:n_series
        if ~mod(j,2) && strcmpi(cfg.graphtype,'cloud_only')
            % move to the other side
            verts{i, j} = [ks{i, j}' + ks_offsets(i) + 0.4, x{i, j}'; ones(n_bins(i, j), 1) * ks_offsets(i) + 0.4, x{i, j}'];
        else
            verts{i, j} = [-ks{i, j}' + ks_offsets(i), x{i, j}'; ones(n_bins(i, j), 1) * ks_offsets(i), x{i, j}'];
        end
    end
end

% fix inputdata for bar plot
limitx = false;
if size(inputdata,2) == 1
    limitx = true;
    inputdata(:,2,1) = 0;
    inputdata(:,2,2) = 0;
end
groupMEANS = squeeze(mean(inputdata,1));

% plot bars
bh = bar(groupMEANS,'EdgeColor','none');
hold on;

pause(.01); % wait to retrieve Xoffset
xoffset(1) = bh(1).XOffset;
xoffset(2) = bh(2).XOffset;

% overwrite plot if irrelevant
if any(strcmpi(cfg.graphtype,{'barcloud', 'bar_only'}))
    for c=1:numel(bh)
        if isfield(cfg,'facecolor')
            if c<=numel(cfg.facecolor)
                set(bh(c),'FaceColor',cfg.facecolor{c});
            end
        end
    end
else
    delete(bh);
end

% plot clouds
if any(strcmpi(cfg.graphtype,{'barcloud','cloud_only'}))
    for i = 1:n_plots_per_series
        for j = 1:n_series
            % plot patches
            h.p{i, j} = patch('Faces', faces{i, j}, 'Vertices', verts{i, j}, 'FaceVertexCData', colours(j, :), 'FaceColor', 'flat', 'EdgeColor', 'none', 'FaceAlpha', cfg.facealpha);
            hold on;
        end
    end
end
hold on;

% offsets for significances
if strcmpi(cfg.graphtype,'barcloud')
    offsettextX = .1;
    offsettextY = cfg.ytickspacing/3;
else
    offsettextX = 0;
    offsettextY = cfg.ytickspacing/2;
end

% h is a cell array containing handles of the various figure parts:
% h.p{i,j} is the handle to the density plot from data{i,j}
% h.s{i,j} is the handle to the 'raindrops' (individual datapoints) from data{i,j}


% plot individual lines
if cfg.plotindiv
    for c = 1:n_plots_per_series
        plot([c+xoffset(1),c+xoffset(2)],squeeze(inputdata(:,c,:)),'-o','color',cfg.linecolor,'MarkerSize',cfg.markersize_indiv);
        hold on;
    end
end

% plot mean dots
if strcmpi(cfg.graphtype, 'cloud_only')
    for c = 1:n_plots_per_series
        plot([c+xoffset(1),c+xoffset(2)],squeeze(mean(inputdata(:,c,:),1)),'-','color','k','LineWidth',2);
        for c2 = 1:n_series
            plot([c+xoffset(c2)],squeeze(mean(inputdata(:,c,c2),1)),'o','color','k','MarkerFaceColor',cfg.facecolor{c2},'MarkerSize',cfg.markersize_mean);
        end
        hold on;
    end
end

if cfg.square
    axis square;
end

if ~isempty(cfg.ylim)
    ylim(cfg.ylim);
    set(gca,'YTick',cfg.ylim(1):cfg.ytickspacing:cfg.ylim(2))
end
xticklabels(cfg.xticklabels);
if ~isempty(cfg.xlabel)
    set(gca,'XTick',1:n_plots_per_series)
    xlabel(cfg.xlabel,'FontSize',cfg.fontsize);
end

if ~isempty(cfg.ylabel)
    ylabel(cfg.ylabel,'FontSize',cfg.fontsize);
end
if limitx
    inputdata(:,2,:) = [];
end
if strcmpi(cfg.graphtype, 'barcloud')
    xlim([0, n_plots_per_series+.5]);
else
    xlim([0.5, n_plots_per_series+.5]);
end

% set tick label sizes
ax = gca;
ax.FontSize = cfg.fontsize; 

% plot p-values
if cfg.plotP
    [~,pVals] = ttest(squeeze(inputdata(:,:,2)),squeeze(inputdata(:,:,1)),'Tail',cfg.tail);
    for c = 1:numel(pVals)
        if pVals(c)<.001
            pCells{c} = 'P<0.001';
        else
            pCells{c} = sprintf('P=%.3f',pVals(c));
        end
    end
    text((1:size(inputdata,2))+offsettextX,max(max(inputdata,[],3))+offsettextY,pCells,'FontSize',cfg.fontsize-2,'horizontalalignment', 'center');
else
    pVals = [];
end
if ~isempty(cfg.conditions)
    leg = legend(cfg.conditions,'Location',cfg.legendlocation,'FontSize',cfg.fontsize);
    leg.ItemTokenSize = [18,18];
    legend boxoff;
end
if ~isempty(cfg.title)
    title(cfg.title,'FontSize',cfg.fontsize+2);
end

% make bar graph translucent and plot on top
if any(strcmpi(cfg.graphtype,{'barcloud', 'bar_only'}))
    for c=1:numel(bh)
        uistack(bh(c),'top');
        bh(c).FaceAlpha = cfg.facealpha;
    end
end

set(gca,'TickDir','out');
set(gca, 'box', 'off');