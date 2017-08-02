function csvstring = cellarray2csvstring(cellarray)
% converts a cell array to a csv string, removes empty values
if size(cellarray,1) == 1
    cellarray = cellarray';
end
cellarray = cellarray(~cellfun(@isempty,cellarray));
cellarray = [cellarray,[repmat({','},numel(cellarray)-1,1);{[]}]]';
csvstring = [cellarray{:}];

