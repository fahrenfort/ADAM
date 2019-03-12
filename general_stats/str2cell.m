function out = str2cell(in,delimiter)
% Cuts a string into a cell array, inserting cell breaks after every delimiter.
% If the delimiter is not present in the string, the function returns the same string.
% J.J.Fahrenfort, VU 2108

if ischar(in)
    breakind = regexp(in,delimiter);
    if ~isempty(breakind)
        breakind = [1 breakind];
        out= cell(1,numel(breakind));
        for cBr = 2:numel(breakind)
            out{cBr} = in(breakind(cBr-1):breakind(cBr));
        end
        out{cBr+1} = in(breakind(end)+1:end);
    end
end
if ~exist('out','var')
    out = in;
end