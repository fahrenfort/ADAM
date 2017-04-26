%% create images for OP experiment
imageDir = '/Users/VU-MBP/Dropbox/Work/Onderwijs/STAGES EN THESES/PSR/capture_masking/pictures/completecars';
outputDir = '/Users/VU-MBP/Dropbox/Work/Onderwijs/STAGES EN THESES/PSR/capture_masking/pictures/newcars';
if ~exist(outputDir,'file')
    mkdir(outputDir);
end

%baseImage = 'Bodies';
baseImage = 'Cars';

picnrs = [118:240];
for c = 1:numel(picnrs)
    % read image and convert to grayscale
    normIm = rgb2gray(imread([imageDir filesep baseImage '_' num2str(picnrs(c)) '.jpg']));
    
    rotIm = rotateImage(normIm,0,400);

    imshow(onlyCircle(rotIm));
    imwrite(onlyCircle(rotIm),[outputDir filesep baseImage '_' num2str(picnrs(c)) '.jpg'],'jpg');
    
    % click to continue
    uiwait(msgbox('Click OK to continue','Waiting...','modal'));
    
%     % scramble phase in 10 steps
%     phaseIm = phaseScramble(rotIm{1},[0:1/19:1],800);
%     for cIm = 1:(numel(phaseIm)-1)
%         imwrite(onlyCircle(phaseIm{cIm}),[outputDir baseImage '_' num2str(c) '_phase' num2str(cIm) '.jpg'],'jpg');
%     end
%     imwrite(onlyCircle(phaseIm{end}),[outputDir baseImage '_' num2str(c) '_randPhase.jpg'],'jpg');
%     
%     % scramble tiles
%     tileSize = 8;
%     tileIm = tileScramble(rotIm{1},tileSize);
%     imwrite(onlyCircle(tileIm),[outputDir baseImage '_' num2str(c) '_randTile' num2str(tileSize) '.jpg'],'jpg');
%     imwrite(tileIm,[outputDir baseImage '_' num2str(c) '_SQrandTile' num2str(tileSize) '.jpg'],'jpg');
%     tileSize = 16;
%     imshow(onlyCircle(tileIm));
%     tileIm = tileScramble(rotIm{1},tileSize);
%     imwrite(onlyCircle(tileIm),[outputDir baseImage '_' num2str(c) '_randTile' num2str(tileSize) '.jpg'],'jpg');
%     imwrite(tileIm,[outputDir baseImage '_' num2str(c) '_SQrandTile' num2str(tileSize) '.jpg'],'jpg');
%     tileSize = 20;
%     imshow(onlyCircle(tileIm));
%     tileIm = tileScramble(rotIm{1},tileSize);
%     imwrite(onlyCircle(tileIm),[outputDir baseImage '_' num2str(c) '_randTile' num2str(tileSize) '.jpg'],'jpg');
%     imwrite(tileIm,[outputDir baseImage '_' num2str(c) '_SQrandTile' num2str(tileSize) '.jpg'],'jpg');
%     tileSize = 25;
%     imshow(onlyCircle(tileIm));
%     tileIm = tileScramble(rotIm{1},tileSize);
%     imwrite(onlyCircle(tileIm),[outputDir baseImage '_' num2str(c) '_randTile' num2str(tileSize) '.jpg'],'jpg');
%     imwrite(tileIm,[outputDir baseImage '_' num2str(c) '_SQrandTile' num2str(tileSize) '.jpg'],'jpg');
%     tileSize = 32;
%     tileIm = tileScramble(rotIm{1},tileSize);
%     imwrite(onlyCircle(tileIm),[outputDir baseImage '_' num2str(c) '_randTile' num2str(tileSize) '.jpg'],'jpg');
%     imwrite(tileIm,[outputDir baseImage '_' num2str(c) '_SQrandTile' num2str(tileSize) '.jpg'],'jpg');
%     tileSize = 40;
%     tileIm = tileScramble(rotIm{1},tileSize);
%     imwrite(onlyCircle(tileIm),[outputDir baseImage '_' num2str(c) '_randTile' num2str(tileSize) '.jpg'],'jpg');
%     imwrite(tileIm,[outputDir baseImage '_' num2str(c) '_SQrandTile' num2str(tileSize) '.jpg'],'jpg');
%     tileSize = 50;
%     tileIm = tileScramble(rotIm{1},tileSize);
%     imwrite(onlyCircle(tileIm),[outputDir baseImage '_' num2str(c) '_randTile' num2str(tileSize) '.jpg'],'jpg');
%     imwrite(tileIm,[outputDir baseImage '_' num2str(c) '_SQrandTile' num2str(tileSize) '.jpg'],'jpg');
% 
%     
    % click to continue
    %uiwait(msgbox(['This was image ' baseImage '_' num2str(c) '. Click OK to move on to the next image.'],'Waiting...','modal'));
    close all;
end