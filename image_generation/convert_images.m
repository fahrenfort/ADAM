%% small script to convert images in the subdirectors of the present
% directory from pgm to png
clear;
dirnames = dir;
dirnames = {dirnames(:).name};
mkdir('png_images');
for cd = 1:numel(dirnames)
    files = dir([pwd filesep dirnames{cd} filesep '*.pgm']);
    files = {files(:).name};
    for cf = 1: numel(files)
        %[pwd filesep dirnames{cd} filesep files{cf}]
        %[pwd filesep 'png_images' filesep files{cf}(1:end-3) 'png']
        imwrite(imread([pwd filesep dirnames{cd} filesep files{cf}]),[pwd filesep 'png_images' filesep files{cf}(1:end-3) 'png']);
    end
end

%% small script to convert jpg to b/w image and to png format
clear;
files = dir('*.jpg');
files = {files(:).name};
for cf = 1: numel(files)
    %[pwd filesep dirnames{cd} filesep files{cf}]
    %[pwd filesep 'png_images' filesep files{cf}(1:end-3) 'png']
    imwrite(rgb2gray(imread([pwd filesep files{cf}])),[pwd filesep files{cf}(1:end-3) 'png']);
end

%% script to make background transparant by creating alpha channel of jpg img and save as png
clear;
files = dir('house*.jpg');
files = {files(:).name};
for cf = 1:numel(files)
    disp(['converting ' pwd filesep files{cf}]);
    image = imread([pwd filesep files{cf}]);
    if size(image,3)>1
        image = rgb2gray(image);
    end
    transval = image(1,1);
    alphamask = ones(size(image));
    alphamask(image>=transval-2 & image<=transval+2) = 0; % how lenient is the algorithm
    L = bwlabel(~alphamask,4);
    labels = unique(L(L>0));
    % only take the backround
    alphamask = ones(size(image));
    for c=1:numel(labels)
        if numel(find(L==labels(c)))>1000 % this indicates the minimum surface to replace
            alphamask(L==labels(c))=0;
        end
    end
    %imshow(alphamask);
    imwrite(image,[pwd filesep files{cf}(1:end-3) 'png'],'png','Alpha',alphamask);
end