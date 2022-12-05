function output = grow_mask(mask, nEl, direction)
% grows 2D mask with non-zero numbers by nEl elements. nEl specifies by how
% many elements the regions should grow by (default: 1).
% direction indicates whether to grow on both sides or just in one
% direction: 'left' (default) 'right' or 'both' (expansion in both
% directions). If direction is 'both' and nEl is an uneven number, the
% algorithm always starts expanding initially towards the left.
% The mask never changes size (overflow is discarded).
% Growing is stopped when a region touches other regions or exceeds the
% dimension of the mask.
% Examples: 
% mask =   [ 0 0 1 1 0 0 0 0 0 1 1 1 0]
% grow_mask(mask, 1)
% becomes: [ 0 1 1 1 0 0 0 0 1 1 1 1 0]
% grow_mask(mask, 3, 'both')
% becomes: [ 0 1 1 1 1 1 0 0 1 1 1 1 1]
%
% by JJF (VU, 2022)

if nargin < 2
    nEl = 1;
end
if nargin < 3
    direction = 'left';
end

rois = bwlabel(mask,8);
nRois = max(unique(rois));
for cRoi=1:nRois
    roiIndex = find(rois==cRoi);
    if strcmpi(direction,'both')
        minIndex = min(roiIndex)-floor(nEl/2);
        maxIndex = max(roiIndex)+ceil(nEl/2);
    elseif strcmpi(direction,'right')
        minIndex = min(roiIndex);
        maxIndex = max(roiIndex)+nEl;
    else
        minIndex = min(roiIndex)-nEl;
        maxIndex = max(roiIndex);
    end
    if minIndex < 1
        minIndex = 1;
    end
    if cRoi > 1 && find(rois==cRoi-1,1,'last') >= minIndex % keep growth from overflowing
        minIndex = find(rois==cRoi-1,1,'last') + 1;
    end
    if maxIndex > numel(mask)
        maxIndex = numel(mask);
    end
    if cRoi < nRois && find(rois==cRoi+1,1,'first') <= maxIndex % keep growth from overflowing
        maxIndex = find(rois==cRoi+1,1,'first') - 1;
    end
    mask(minIndex:maxIndex) = 1;
end
output = mask;
