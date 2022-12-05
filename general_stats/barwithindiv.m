function pVals = barwithindiv(cfg,inputdata)
% function barwithindiv(cfg,inputdata)
% Makes a grouped bar plot with individual inputdata points connected to
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
% defcfg.xlabel = '';               % x-axis label
% defcfg.ylabel = '';               % y-axis label
% defcfg.title = '';                % title of plot
% defcfg.plotP = true;              % perform t-test and plot p-values, p-values are also returned by the function 
% cfg.tail = 'both';                % indicate that t-test should be double sided, alternatively: 'right' or 'left' (see help ttest) 
% defcfg.square = true;             % plot the axes as square or not


% defaults
defcfg.linecolor = [.4 .4 .4];
defcfg.conditions = [];
defcfg.legendlocation = 'NorthEast';
defcfg.facecolor{1} = [0 118 186]/255; % colors of the bars
defcfg.facecolor{2} = [181 23 0]/255;
defcfg.xlim = [];
defcfg.ylim = [];
defcfg.xticklabels = [];
defcfg.xlabel = '';
defcfg.ylabel = '';
defcfg.title = '';
defcfg.plotP = true;    % perform t-test and plot p-values
defcfg.tail = 'both';   % whether t-test is two-sided or one-sided (see help ttest) 
defcfg.square = false;   % plot the axes as square

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

% fix inputdata
limitx = false;
if size(inputdata,2) == 1
    limitx = true;
    inputdata(:,2,1) = 0;
    inputdata(:,2,2) = 0;
end

groupMEANS = squeeze(mean(inputdata,1));
% groupSE = squeeze(std(inputdata,1)/sqrt(size(inputdata,1)));
bh = bar(groupMEANS);
%bh = barwitherr(groupSE,groupMEANS);
if cfg.square
    axis square;
end
% plot individual lines
for c = 1:size(inputdata,2)
    hold on;
    plot([c+bh(1).XOffset c+bh(2).XOffset],squeeze(inputdata(:,c,:)),'-o','color',cfg.linecolor);
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
end
xticklabels(cfg.xticklabels);
if ~isempty(cfg.xlabel)
    xlabel(cfg.xlabel);
end
if ~isempty(cfg.ylabel)
    ylabel(cfg.ylabel);
end
if limitx
    inputdata(:,2,:) = [];
    xlim([.5, 1.5]);
end
% plot p-values
if cfg.plotP
    [~,pVals] = ttest(squeeze(inputdata(:,:,2)),squeeze(inputdata(:,:,1)),'Tail',cfg.tail);
    text((1:size(inputdata,2))-.1,max(max(inputdata,[],3))+.02,num2cell(round(pVals(:),3)));
else
    pVals = [];
end
if ~isempty(cfg.conditions)
    legend(cfg.conditions,'Location',cfg.legendlocation);
    legend boxoff;
end
if ~isempty(cfg.title)
    title(cfg.title,'FontSize',12);
end
for c=1:numel(bh)
    uistack(bh(c),'top');
    bh(c).FaceAlpha = 0.7;
end
