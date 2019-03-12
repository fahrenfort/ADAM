function out = str2cell(in,delimiter)
% function out = str2cell(in,delimiter)
% Cuts a string into a cell array, inserting cell breaks at every delimiter. The default delimiter
% is a comma (','). If the delimiter is not present in the string, the function returns the same
% input string without converting to cell.
% J.J.Fahrenfort, VU 2018, 2019
if nargin<2
    delimiter = ',';
end
if ischar(in)
    breakind = regexp(in,delimiter,'once');
    if ~isempty(breakind)
        if ~strcmp(in(1),delimiter)
            in = [delimiter in];
        end
        if ~strcmp(in(end),delimiter)
            in = [ in delimiter ];
        end
        breakind = regexp(in,delimiter);
        out = cell(1,numel(breakind)-1);
        for cBr = 1:numel(breakind)-1
            out{cBr} = in(breakind(cBr)+1:breakind(cBr+1)-1);
        end
    end
end
if ~exist('out','var')
    out = in;
end