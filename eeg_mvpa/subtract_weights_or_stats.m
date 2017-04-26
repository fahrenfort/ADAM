function dif_out = subtract_weights_or_stats(statsorweights1,statsorweights2)
% function difstatsorweights = subtract_MVPA_conditions(statsorweights1,statsorweights2)
% subtracts weight or stats structs from each other according to
% difstatsorweights = statsorweights1 - statsorweights2
% 
structfields = fieldnames(statsorweights1);
for c = 1:numel(structfields)
    if sum(strcmpi(structfields{c},{'chanlocs','channelpool','settings'}))>0
        dif_out.(structfields{c}) = statsorweights1.(structfields{c});
    elseif any(strcmpi(structfields{c},{'condname'}))
        dif_out.(structfields{c}) = [ statsorweights1.condname '-' statsorweights2.condname ];
    elseif any(strcmpi(structfields{c},{'CTFpercond','semCTFpercond','indivCTFpercond'}))
        % we'll get to this later
    elseif strcmpi(structfields{c},'indivCTF')
        dif_out.(structfields{c}) = statsorweights1.(structfields{c}) - statsorweights2.(structfields{c});
        dif_out.semCTF = std(dif_out.(structfields{c}))/sqrt(size(dif_out.(structfields{c}),1));
    elseif any(strcmpi(structfields{c},{'StdError','pVals','pStruct'}))
        % these need to be recomputed by the plot function operating on the
        % output of the current function
    elseif ~iscell(statsorweights1.(structfields{c})) && ~strcmpi(structfields{c},'semCTF')
        dif_out.(structfields{c}) = statsorweights1.(structfields{c}) - statsorweights2.(structfields{c});
    end
end
