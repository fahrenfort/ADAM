function nt_topoplot(cfg,data)
%nt_topoplot(cfg,data) - simple topoplot
%
%  nt_topoplot(topo,data): use 'topo' layout file
%

if nargin<2; error('!'); end

data=data(:);

if ~isstruct(cfg)
    topo=cfg; 
    cfg=[];
    cfg.layout=topo;
end

% todo: extract fields from cfg as arguments to ft_topo_plot
       
layout=ft_prepare_layout(cfg);


ft_plot_topo(layout.pos(1:numel(data),1), layout.pos(1:numel(data),2), data',...
    'mask', layout.mask,...
    'outline', layout.outline, ...
    'interplim', 'mask');

% ensure symmetrical color axes
clim=get(gca,'clim');
set(gca,'clim',[-1 1]*max(abs(clim)));

axis off
axis equal
