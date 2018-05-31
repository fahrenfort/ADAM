function [ colorblind_matlab, colorblind_palette ] = color_blindness_palette()
% make color palette for the color blind
% useful for figures that are readable for everybody
% see color_blindness_palette.png in this folder
colorblind_palette = {[0 0 0],    [0 73 73],    [0 146 146],   [255 109 182], [255 182 219] ...
                       [73 0 146], [0 109 219], [182 109 255], [109 182 255], [182 219 255] ...
                       [146 0 0],  [146 73 0],  [219 109 0],   [36 255 36],   [255 255 109]};
colorblind_matlab = cellfun(@(x) x/255, colorblind_palette, 'UniformOutput',false);