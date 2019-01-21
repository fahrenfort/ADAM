function nt_mark(idx,labels,line_params,text_params)
% nt_mark(idx,labels,line_params,text_params)


if nargin<2; labels=[]; end

ylim=get(gca,'ylim');
for k=1:numel(idx)
    line([idx(k) idx(k)],[ylim(1),ylim(1)+0.95*(ylim(2)-ylim(1))], 'color', 'r', 'linestyle', '--');
    if ~isempty(labels)
        if k<= numel(labels) && ~isempty(labels{k});
            text(idx(k), ylim(1)+0.97*(ylim(2)-ylim(1)), labels{k});
        end
    end
end
set(gca,'ylim',ylim)
