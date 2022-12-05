function output = und2spc(input)
% replace underscores with spaces
if iscell(input)
    output = cell(size(input));
    for c = 1:numel(input)
        output{c} = und2spc(input{c});
    end
else
    output = strrep(input,'_',' ');
end