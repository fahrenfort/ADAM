function csvstring = cell2csv(cellarray)
% Convert cell array into a single csv string (columnwise). Returns input if it's already a string.
% If input array is a double, converts it to a cell array. Removes empty values, converts numerics
% to chars. Useful when wanting to concisely output text to the screen using disp() or otherwise.
%
% Example usage:
%
% cell2csv({'aap' 'noot' 'mies' 1 2 3})
%
% ans =
% aap,noot,mies,1,2,3
%
% By J.J.Fahrenfort, VU/UvA 2018

if ischar(cellarray)
    csvstring = cellarray;
    return
end
cellarray = shiftdim(cellarray);
if isnumeric(cellarray)
    cellarray = num2cell(cellarray);
end
cellarray = cellfun(@num2str,cellarray,'UniformOutput',false);
cellarray = cellarray(~cellfun(@isempty,cellarray));
cellarray = [cellarray,[repmat({','},numel(cellarray)-1,1);{[]}]]';
csvstring = [cellarray{:}];

