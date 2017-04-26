function new_weights = subtract_FEM_channels(femweights,condition_def)
% function new_weights = subtract_FEM_channels(femweights,condition_def)
% subtract weights in channel A from weights in channel B, as outlined in
% condition_def{1} = [B,A];
% condition_def is a cell array of condition definitions such that
% condition_def{1} = [2,3]; 
% condition_def{2} = [4,5]; 
% means that a new set of FEM weights will be generated in which condition
% 3 was subtracted from condition 2 and a second set of weights in which
% condition 5 was subtracted from condition 4.
% output weights can be plotted using plot_FEM_weights

for c = 1:numel(femweights)
    new_weights(c) = sub_subtract_weigths(femweights(c),condition_def);
end

function new_weights = sub_subtract_weigths(weights,condition_def)
for c = 1:numel(condition_def)
    indivWeights(:,:,:,c) = weights.indivWeights(:,:,:,condition_def{c}(1)) - weights.indivWeights(:,:,:,condition_def{c}(2));
    avWeights(:,:,c) = weights.avWeights(:,:,condition_def{c}(1)) - weights.avWeights(:,:,condition_def{c}(2));
end
new_weights.indivWeights = indivWeights;
new_weights.avWeights = avWeights;
new_weights.chanlocs = weights.chanlocs;