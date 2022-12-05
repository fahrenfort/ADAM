function output = spc2und(input)
% replace underscores with spaces
if iscell(input)
    output = cell(size(input));
    for c = 1:numel(input)
        output{c} = spc2und(input{c});
    end
else
    output = strrep(input,' ','_');
end