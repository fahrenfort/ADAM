% This script shows how p-values and Bayes factors behave under sequential
% testing when there is no effect in your data. It requires the Bayes
% Factor package: https://www.mathworks.com/matlabcentral/fileexchange/69794-bayesfactor
cFig = 0;
figure; 
for it = 1:10
    disp(['simulation ' num2str(it)]); % progress bar
    cStat = 0;
    clear BFs pVals;
    nullData = randn([2,100000]);
    for c=1:100:size(nullData,2)
        cStat = cStat + 1;
        % running sequentially (on cumuluative data):
        [BFs(cStat), pVals(cStat)] = bf.ttest(nullData(1,1:c),nullData(2,1:c)); 
        % running on new samples:
        % [BFs(cStat), pVals(cStat)] = bf.ttest(nullData(1,c:c+98),nullData(2,c:c+98));
    end
    % 1
    subplot(5,8,1+cFig);
    hold on;
    yline(.05,'--');
    plot(pVals);
    ylim([0,1]);
    ntitle('pvals');
    % 2
    subplot(5,8,2+cFig);
    histogram(pVals,'BinWidth',0.05);
    xlim([0,1]);
    ntitle('hist pvals');
    FAR = sum(pVals<.05)/numel(pVals); % false alarm rate (proportion of P<.05)
    text(max(get(gca,'XLim'))/2,3*max(get(gca,'YLim'))/4,sprintf('FAR=%1.2f',FAR));
    % 3
    subplot(5,8,3+cFig);
    plot(BFs,'r');
    hold on;
    yline(1,'--');
    ylim([0,3]);
    ntitle('BFs');
    % 4
    subplot(5,8,4+cFig);
    histogram(BFs,'FaceColor','r','BinWidth',1/6);
    xlim([0,3]);
    ntitle('hist BFs');
    BF0 = sum(BFs<.33)/numel(BFs); % proportion of BFs for the null (BF<.33)
    BF1 = sum(BFs>3)/numel(BFs); % proportion of BFs for H1 (BF>3)
    text(max(get(gca,'XLim'))/8,3*max(get(gca,'YLim'))/4,sprintf('BF0=%1.2f/BF1=%1.2f',BF0,BF1));
    %
    cFig = cFig + 4;
end

