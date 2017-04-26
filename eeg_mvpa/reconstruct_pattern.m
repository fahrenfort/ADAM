function pattern = reconstruct_pattern(channel_response, femweights)
% reconstruct a new pattern based on a channel response

% run through stats
for c = 1:numel(femweights)
    pattern(c) = sub_reconstruct_pattern(channel_response, femweights(c));
end

function pattern = sub_reconstruct_pattern(channel_response, femweights)
for cSubj = 1:size(femweights.indivWeights,1)
    for cT = 1:size(femweights.indivWeights,2)
        W = squeeze(femweights.indivWeights(cSubj,cT,:,:));
        indivPattern(cSubj,cT,:,:) = W*channel_response;
    end
end
avPattern = squeeze(mean(indivPattern,1));
pattern.indivWeights = indivPattern;
pattern.avWeights = avPattern;
pattern.chanlocs = femweights.chanlocs;
