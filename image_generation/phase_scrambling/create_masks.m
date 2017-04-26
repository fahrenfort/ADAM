%% create attneave shape for mask
source_dir = '/Users/VU-MBP/Dropbox/Work/experimenten/Attention_masked/capture_masking/pictures/masks';
target_dir = '/Users/VU-MBP/Dropbox/Work/experimenten/Attention_masked/new_masks';

rand_series = randperm(240);
nColors = 128; % max 128
for c = 1:2:240
    [mask1, map1] = rgb2ind(imread([source_dir filesep 'Masks_' num2str(rand_series(c)) '.jpg'],'jpg'),nColors);
    [mask2, map2] = rgb2ind(imread([source_dir filesep 'Masks_' num2str(rand_series(c+1)) '.jpg'],'jpg'),nColors);
    contrast = rand(1)/5+.5; % scales it between .5 and .7
    map1 = rgb2gray(brighten(map1,+contrast));
    map2 = rgb2gray(brighten(map2,-contrast));
    ImageSize = [400 400];
    mask1 = mask1(1:ImageSize(1),1:ImageSize(2));
    mask2 = mask2(1:ImageSize(1),1:ImageSize(2));
    [~, ~, Xp, Yp]=ShapeFamily('AngleLims',[20 175],'NSides', 20, 'NPts2Shift', 20, 'NMembers', 1,'MakePix','n');
    mask_image = uint8(zeros(ImageSize(1), ImageSize(2)));
    im_bool = flipud(~roipoly(mask_image, Xp * ImageSize(2), Yp * ImageSize(1)));
    mask_image(im_bool) = uint8(mask1(im_bool));
    mask_image(~im_bool) = uint8(mask2(~im_bool))+nColors;
    imwrite(mask_image,[map1;map2],[target_dir filesep 'ObjMasks_' num2str(c) '.jpg'],'jpg');
    mask_image(~im_bool) = uint8(mask1(~im_bool));
    mask_image(im_bool) = uint8(mask2(im_bool))+nColors;
    imwrite(mask_image,[map1;map2],[target_dir filesep 'ObjMasks_' num2str(c+1) '.jpg'],'jpg');
    %imshow(mask_image,[map1;map2]);
end
