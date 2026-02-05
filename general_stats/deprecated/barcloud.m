function pVals = cloudwithindiv(cfg,inputdata)
% function pVals = cloudwithindiv(cfg,inputdata)
% Makes a grouped raincloud plot with individual inputdata points connected to
% each other, including p-values above each bar combination. Can be useful
% when plotting data from a 2-factor factorial design.
% inputdata must have the following format:
% inputdata(participants, rep_meas1, rep_meas2);
% rep_meas2 must have at most 2 levels if one wants to perform a t-test
% between those levels
% Various optional arguments can be given to adjust colors, axis names, xtick labels
% using the form:
% defcfg.linecolor = [.4 .4 .4];    % to indicate the color of the lines of individual subjects
% defcfg.conditions = {'1', '2'};   % indicate condition names that the bars refer to
% cfg.facecolor{1} = [0 0.46 0.73]; % colors of the bars
% cfg.facecolor{2} = [0.73 0.46 0];
% defcfg.legendlocation = 'NorthEast'; % where to put the legend (see help legend) 
% defcfg.xlim = [];                 % limits on x-axis
% defcfg.ylim = [];                 % limits on y-axis
% defcfg.xticklabels = [];          % labels on levels for conditions
% defcfg.ytickspacing = .1;         % spacing between ticks on the y-axis
% defcfg.xlabel = '';               % x-axis label
% defcfg.ylabel = '';               % y-axis label
% defcfg.fontsize = 12;
% defcfg.title = '';                % title of plot
% defcfg.plotP = true;              % perform t-test and plot p-values, p-values are also returned by the function 
% defcfg.tail = 'both';                % indicate that t-test should be double sided, alternatively: 'right' or 'left' (see help ttest) 
% defcfg.square = true;                % plot the axes as square or not
% defcfg.FaceAlpha = .55;
%
% J.J.Fahrenfort, VU 2022

% defaults
defcfg.linecolor = [.4 .4 .4];
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
defcfg.FaceAlpha = .55;
defcfg.density_granularity = 200;

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

% Calculate properties of density plots
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

% determine spacing between plots
%spacing     = 2 * mean(mean(cellfun(@max, ks)));
normalize = num2cell(repmat(2.2 * max(max(cellfun(@max, ks))),size(ks)));
ks = cellfun(@mrdivide,ks,normalize,'UniformOutput',false);

% if max(max(cellfun(@max, ks)))>.45
%     disp('help');
% end

spacing = 1;
ks_offsets  = [.7:n_plots_per_series] .* spacing;

% calculate patch vertices from kernel density
for i = 1:n_plots_per_series
    for j = 1:n_series
        verts{i, j} = [-ks{i, j}' + ks_offsets(i), x{i, j}'; ones(n_bins(i, j), 1) * ks_offsets(i), x{i, j}'];
        verts{i, j} = [-ks{i, j}' + ks_offsets(i), x{i, j}'; ones(n_bins(i, j), 1) * ks_offsets(i), x{i, j}'];        
    end
end

% patches
for i = 1:n_plots_per_series
    for j = 1:n_series
        % plot patches
        h.p{i, j} = patch('Faces', faces{i, j}, 'Vertices', verts{i, j}, 'FaceVertexCData', colours(j, :), 'FaceColor', 'flat', 'EdgeColor', 'none', 'FaceAlpha', cfg.FaceAlpha);
    end
end
hold on;

% fix inputdata
limitx = false;
if size(inputdata,2) == 1
    limitx = true;
    inputdata(:,2,1) = 0;
    inputdata(:,2,2) = 0;
end

% bh = rm_raincloud(cloud_data,facecolor,0);

groupMEANS = squeeze(mean(inputdata,1));

% plot bars
bh = bar(groupMEANS,'EdgeColor','none');

% h is a cell array containing handles of the various figure parts:
% h.p{i,j} is the handle to the density plot from data{i,j}
% h.s{i,j} is the handle to the 'raindrops' (individual datapoints) from data{i,j}

pause(.01); 

% plot individual lines
if cfg.plotindiv
    for c = 1:n_plots_per_series
        hold on;
        plot([c+bh(1).XOffset c+bh(2).XOffset],squeeze(inputdata(:,c,:)),'-o','color',cfg.linecolor);
    end
end
hold on;

if cfg.square
    axis square;
end

for c=1:numel(bh)
    if isfield(cfg,'facecolor')
        if c<=numel(cfg.facecolor)
            set(bh(c),'FaceColor',cfg.facecolor{c});
        end
    end
end
if ~isempty(cfg.xlim)
    xlim(cfg.xlim);
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
xlim([0, n_plots_per_series+.5]);

% set tick label sizes
ax = gca;
ax.FontSize = cfg.fontsize; 

% plot p-values
if cfg.plotP
    [~,pVals] = ttest(squeeze(inputdata(:,:,2)),squeeze(inputdata(:,:,1)),'Tail',cfg.tail);
    for c = 1:numel(pVals)
        pCells{c} = sprintf('P=%.3f',pVals(c));
    end
    text((1:size(inputdata,2))-.2,max(max(inputdata,[],3))+.02,pCells,'FontSize',cfg.fontsize-2);
    %text(max(max(inputdata,[],3))+.02,(1:size(inputdata,2))-.1,num2cell(round(pVals(:),3)));
else
    pVals = [];
end
if ~isempty(cfg.conditions)
    leg = legend(cfg.conditions,'Location',cfg.legendlocation,'FontSize',cfg.fontsize);
    leg.ItemTokenSize = [18,18];
    legend boxoff;
end
if ~isempty(cfg.title)
    title(cfg.title,'FontSize',cfg.fontsize);
end

% lower point size
% h.s{i,j} is the handle to the 'raindrops' (individual datapoints) from data{i,j}
% h.m(i,j) is the handle to the single, large dot that represents mean(data{i,j})
% h.l(i,j) is the handle for the line connecting h.m(i,j) and h.m(i+1,j)
for c1 = 1:size(inputdata,2)
    for c2 = 1:size(inputdata,3)
%        bh.s{c1,c2}.SizeData = 50;
%        bh.m(c1,c2).SizeData = 75;
    end
end

for c=1:numel(bh)
    uistack(bh(c),'top');
    bh(c).FaceAlpha = cfg.FaceAlpha;
end

set(gca,'TickDir','out');