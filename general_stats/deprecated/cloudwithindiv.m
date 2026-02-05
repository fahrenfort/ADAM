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
% defcfg.xlabel = '';               % x-axis label
% defcfg.ylabel = '';               % y-axis label
% defcfg.title = '';                % title of plot
% defcfg.plotP = true;              % perform t-test and plot p-values, p-values are also returned by the function 
% cfg.tail = 'both';                % indicate that t-test should be double sided, alternatively: 'right' or 'left' (see help ttest) 
% defcfg.square = true;             % plot the axes as square or not
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
defcfg.xlabel = '';
defcfg.ylabel = '';
defcfg.title = '';
defcfg.plotP = true;    % perform t-test and plot p-values
defcfg.tail = 'both';   % whether t-test is two-sided or one-sided (see help ttest) 
defcfg.square = false;   % plot the axes as square
defcfg.plotindiv = false;

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
    facecolor(c,:) = cfg.facecolor{c};
end
for c1 = 1:size(inputdata,2)
    for c2 = 1:size(inputdata,3)
        cloud_data{c1,c2} = squeeze(inputdata(:,c1,c2));
    end
end

% % fix inputdata
% limitx = false;
% if size(inputdata,2) == 1
%     limitx = true;
%     inputdata(:,2,1) = 0;
%     inputdata(:,2,2) = 0;
% end

bh = rm_raincloud(cloud_data,facecolor,0);

% groupMEANS = squeeze(mean(inputdata,1));

% plot bars
% bh = bar(groupMEANS);
%bh = barwitherr(groupSE,groupMEANS);


% h is a cell array containing handles of the various figure parts:
% h.p{i,j} is the handle to the density plot from data{i,j}
% h.s{i,j} is the handle to the 'raindrops' (individual datapoints) from data{i,j}


% % plot individual lines
% if cfg.plotindiv
%     for c = 1:size(inputdata,2)
%         hold on;
%         plot([c+bh(1).XOffset c+bh(2).XOffset],squeeze(inputdata(:,c,:)),'-o','color',cfg.linecolor);
%     end
% end
% for c=1:numel(bh)
%     if isfield(cfg,'facecolor')
%         if c<=numel(cfg.facecolor)
%             set(bh(c),'FaceColor',cfg.facecolor{c});
%         end
%     end
% end
% if ~isempty(cfg.xlim)
%     ylim(cfg.xlim);
% end
if ~isempty(cfg.ylim)
    xlim(cfg.ylim);
end
yticklabels(cfg.xticklabels(end:-1:1));
if ~isempty(cfg.xlabel)
    ylabel(cfg.xlabel);
end
if ~isempty(cfg.ylabel)
    xlabel(cfg.ylabel);
end
% if limitx
%     inputdata(:,2,:) = [];
%     xlim([.5, 1.5]);
% end
% plot p-values
if cfg.plotP
    [~,pVals] = ttest(squeeze(inputdata(:,:,2)),squeeze(inputdata(:,:,1)),'Tail',cfg.tail);
    %text((1:size(inputdata,2))-.1,max(max(inputdata,[],3))+.02,num2cell(round(pVals(:),3)));
    text(max(max(inputdata,[],3))+.02,(1:size(inputdata,2))-.1,num2cell(round(pVals(:),3)));
else
    pVals = [];
end
if ~isempty(cfg.conditions)
    legend([ bh.p{1,1} bh.p{1,2}],cfg.conditions,'Location',cfg.legendlocation);
    %legend([ bh.p{1,1} bh.p{1,2}],cfg.conditions);
    legend boxoff;
end
if ~isempty(cfg.title)
    title(cfg.title,'FontSize',12);
end

% lower point size
% h.s{i,j} is the handle to the 'raindrops' (individual datapoints) from data{i,j}
% h.m(i,j) is the handle to the single, large dot that represents mean(data{i,j})
% h.l(i,j) is the handle for the line connecting h.m(i,j) and h.m(i+1,j)
for c1 = 1:size(inputdata,2)
    for c2 = 1:size(inputdata,3)
        bh.s{c1,c2}.SizeData = 50;
        bh.m(c1,c2).SizeData = 75;
    end
end

if cfg.square
    axis square;
end

% for c=1:numel(bh)
%     uistack(bh(c),'top');
%     bh(c).FaceAlpha = 0.7;
% end
