function nt_same_clim(h)
%nt_same_clim(h) - harmonize color limits of plots within figure
%
%  h: figure handle [default: current]

if nargin<1; h=gcf; end

c=get(h,'children');
clims=[];
for iC=1:numel(c)    
    cc=c(iC);
    try
        if ~isempty(get(cc,'clim'));
            clims=[clims; (get(cc,'clim'))];
        end
    end
end
for iC=1:numel(c)    
    cc=c(iC);
    try
        set(cc,'clim',[min(clims(:)),max(clims(:))]);
    end
end

        
        