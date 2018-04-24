function dif_out = adam_subtract_stats(statsorweights1,statsorweights2)
% function difstatsorweights = subtract_MVPA_conditions(statsorweights1,statsorweights2)
% subtracts weight or stats structs from each other according to
% difstatsorweights = statsorweights1 - statsorweights2
% 
structfields = fieldnames(statsorweights1);
for c = 1:numel(structfields)
    if any(strcmpi(structfields{c},{'chanlocs','channelpool','settings','cfg'}))
        dif_out.(structfields{c}) = statsorweights1.(structfields{c});
    elseif any(strcmpi(structfields{c},{'weights'}))
        dif_out.(structfields{c}) = adam_subtract_stats(statsorweights1.(structfields{c}),statsorweights2.(structfields{c}));
    elseif any(strcmpi(structfields{c},{'condname'}))
        dif_out.(structfields{c}) = [ statsorweights1.condname '-' statsorweights2.condname ];
    elseif any(strcmpi(structfields{c},{'CTFpercond','semCTFpercond','indivCTFpercond'}))
        % we'll get to this later
    elseif strcmpi(structfields{c},'indivCTF')
        dif_out.(structfields{c}) = statsorweights1.(structfields{c}) - statsorweights2.(structfields{c});
        dif_out.semCTF = std(dif_out.(structfields{c}))/sqrt(size(dif_out.(structfields{c}),1));
    elseif any(strcmpi(structfields{c},{'StdError','pVals','pStruct','mpcompcor_method'}))
        dif_out.(structfields{c}) = [];
        % these need to be recomputed by the plot function operating on the
        % output of the current function
    elseif ~iscell(statsorweights1.(structfields{c})) && ~strcmpi(structfields{c},'semCTF')
        disp(['subtracting ' structfields{c}]);
        dif_out.(structfields{c}) = statsorweights1.(structfields{c}) - statsorweights2.(structfields{c});
    end
end
