%% create images for OP experiment
imageDir = '/Users/VU-MBP/Dropbox/Work/Experimenten/Attention_masked/capture_exp/pictures/scenes';
outputDir = '/Users/VU-MBP/Dropbox/Work/Experimenten/Attention_masked/capture_exp/pictures/newprimes';
if ~exist(outputDir,'file')
    mkdir(outputDir);
end

baseImage = 'Bodies';
for c = 1:1
    % read image and convert to grayscale
    normIm = rgb2gray(imread([imageDir filesep baseImage '_' num2str(c) '.jpg']));
    
    % scramble phase, add noise 
    phaseIm = phaseScramble(normIm,[0:1/8:1],800);
    for cIm = 1:numel(phaseIm)        
        imwrite(phaseIm{cIm},[outputDir filesep baseImage '_' num2str(c) '_phase' num2str(cIm) '.jpg'],'jpg');
    end
    
    for cIm = 1:numel(phaseIm)
        %imwrite(imnoise(normIm,'gaussian',0,(log(cIm)*(cIm/numel(phaseIm)))^2),[outputDir baseImage '_' num2str(c) '_noise' num2str(cIm) '.jpg'],'jpg');
        imwrite(imnoise(normIm,'gaussian',0,log(cIm)^4),[outputDir baseImage '_' num2str(c) '_noise' num2str(cIm) '.jpg'],'jpg');
    end

    temp = normIm;
    normIm = double(normIm/255);
    minIntens = double(min(min(normIm)));
    maxIntens = double(max(max(normIm)));
    meanIntens = double(mean(mean(normIm)));
    lowerLimit = [minIntens:(meanIntens-minIntens)/8:meanIntens];
    upperLimit = [maxIntens:-(maxIntens-meanIntens)/8:meanIntens];
    normIm = temp;
    for cIm = 1:numel(lowerLimit)
        imwrite(imadjust(normIm,[],[lowerLimit(cIm) upperLimit(cIm)]),[outputDir baseImage '_' num2str(c) '_contrast' num2str(cIm) '.jpg'],'jpg');
    end
    
%     tileSize = [1 8 16 20 25 32 40 50 100 200 400];
%     for cIm = 1:numel(tileSize)
%         tileIm = tileScramble(normIm,tileSize(cIm));
%         imwrite(tileIm,[outputDir baseImage '_' num2str(c) '_SQrandTile' num2str(tileSize) '.jpg'],'jpg');
%     end
end