function [connectivity] = get_connected_electrodes(chanlocs)
% function to obtain connectivity between electrodes
% uses some Fieldtrip functions (thanks guys)
% output describes for each electrode (rows) to which electrodes it is
% connected (columns)
% JJF, VU 2016

% code below is only needed when the channel positions are not available, which should be fixed at
% the first level analysis and/or when reading in the data rather than here (so this can be deleted,
% just keeping it in for now)
% elec = ft_read_sens('plotting_1005.sfp');
% [~, ~, tokeep] = intersect(labels,elec.label,'stable');
% elec.chanpos = elec.chanpos(tokeep,:);
% elec.elecpos = elec.elecpos(tokeep,:);
% elec.label = elec.label(tokeep);
% elec.type = 'eeg1020';
% if sum(strcmpi(elec.label', labels)) ~= numel(labels)
%     error('hmm, this should not happen. need to fix label order, get to work');
% end

% restructure positions to be able to work with fieldtrip
elec.chanpos = [ chanlocs(:).X ; chanlocs(:).Y ; chanlocs(:).Z]';
elec.elecpos = elec.chanpos;
elec.label = {chanlocs(:).labels};

cfg = [];
cfg.elec = elec;
cfg.method = 'triangulation';
neighbours = ft_prepare_neighbours(cfg);
cfg = [];
cfg.channel = elec.label;
cfg.neighbours = neighbours;
cfg.connectivity = channelconnectivity(cfg);
connectivity = cfg.connectivity;
