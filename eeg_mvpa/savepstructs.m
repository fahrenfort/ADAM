function savepstructs(stats,fname)
% saves all the p-structs in stats as txt file
% J.J.Fahrenfort, VU 2017
if isfield(stats(1),'pStruct') && ~isempty([stats(:).pStruct])
    clustertypes = fieldnames(stats(1).pStruct);
    fid = fopen([fname '.txt'],'w');
    for cStats = 1:numel(stats)
        fprintf(fid,'%s\n\n',stats(cStats).condname);
        for cType = 1:numel(clustertypes)
            clusters = stats(cStats).pStruct.(clustertypes{cType});
            for cClust =1:numel(clusters)
                flds = fieldnames(clusters(cClust));
                for cFld = 1:numel(flds)
                    fprintf(fid,'%s:\t%g\n',flds{cFld},clusters(cClust).(flds{cFld}));
                end
                fprintf(fid,'\n');
            end
            fprintf(fid,'\n');
        end
    end
else
    disp('no pstructs to save to file');
end