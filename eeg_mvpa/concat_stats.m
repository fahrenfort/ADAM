function outstats = concat_stats(varargin)
% ADAM_CONCAT_STATS concatenates stats structures if they do not have the same fields by removing
% fields that do not exist in both structures.
%
% Usage: outstats = adam_concat_stats(stats1, stats2, ...)
% By J.J.Fahrenfort, VU, 2018

allflds = fieldnames(varargin{1});
for cStats = 2:numel(varargin)
    flds = fieldnames(varargin{cStats});
    allflds = intersect(allflds, flds);
end

outstats = [];
for cStats = 1:numel(varargin)
    flds = fieldnames(varargin{cStats});
    tmp = varargin{cStats};
    for cFlds = 1:numel(flds)
        if ~ismember(flds{cFlds},allflds)
            tmp = rmfield(tmp,flds{cFlds});
        end
    end
    outstats = [outstats tmp];
end