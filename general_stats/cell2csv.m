function csvstring = cell2csv(cellarray,addspace)
% Convert cell array into a single comma separated value string (csv). Operates columnwise. Returns
% input if it's already a string. Removes empty values, converts numerics to chars. Useful when
% wanting to concisely output text to the screen using disp() or otherwise.
%
% Inputs:
%           cellarray   - cell array of strings. If cellarray is a double, converts it to a cell 
%                         array of strings prior to creating the csv. 
%           addspace    - false (default) or true. If true, adds a space to every comma in the csv.
% Example usage:
%
% cell2csv({'aap' 'noot' 'mies' 1 2 3})
%
% ans =
% aap,noot,mies,1,2,3
%
% By J.J.Fahrenfort, VU/UvA 2018

if nargin < 2
    addspace = false;
end

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
if addspace
    cellarray = [cellarray,[repmat({', '},numel(cellarray)-1,1);{[]}]]';
else
    cellarray = [cellarray,[repmat({','},numel(cellarray)-1,1);{[]}]]';
end
csvstring = [cellarray{:}];

