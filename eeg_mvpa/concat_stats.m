function outstats = concat_stats(varargin)
% ADAM_CONCAT_STATS concatenates stats structures if they do not have the same fields by removing
% fields that do not exist in both structures.
%
% Usage: outstats = concat_stats(stats1, stats2, ...)
% By J.J.Fahrenfort, VU, 2018

allflds = fieldnames(varargin{1});
for cStats = 2:numel(varargin)
    flds = fieldnames(varargin{cStats});
    allflds = intersect(allflds, flds);
end

% remove empty fields
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

% re-order to have difference stats last
if isfield(outstats(1),'settings')
    moveToEnd = [];
    moveToStart = [];
    for cStats = 1:numel(outstats)
        if ~isempty(strfind(outstats(cStats).settings.measuremethod,' difference')) || ~isempty(strfind(outstats(cStats).settings.measuremethod,' correlation'))
            moveToEnd = [moveToEnd cStats];
        else
            moveToStart = [moveToStart cStats];
        end
    end
    outstats = [outstats(moveToStart) outstats(moveToEnd)];
end