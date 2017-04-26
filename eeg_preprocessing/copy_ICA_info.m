function copy_ICA_info(filepath,filenames_in,filenames_out)
% function copy_ICA_info(filepath,filenames_in,filenames_out)
% copies all the ICA information present in filenames_in to filenames_out
% run from within matlab
% filepath: where files are stored
% filenames_in: cell array with names of files
% filenames_out: cell array with names of files
% filenames_in and filenames_out should contain the same number of files
% Example usage:
% copy_ICA_info(pwd,{'subj1_task1' 'subj2_task1'},{'subj1_task2' 'subj2_task2'});

if numel(filenames_in) ~= numel(filenames_out)
    disp('number of files in filenames_in and filenames_out does not correspond')
    help copy_ICA_info;
    return;
end

for cFiles = 1:numel(filenames_in)
    filename_in = filenames_in{cFiles};
    filename_out = filenames_out{cFiles};
    if ~strcmp(filename_in(end-2:end),'set')
        filename_in = [filename_in '.set'];
    end
    if ~strcmp(filename_out(end-2:end),'set')
        filename_out = [filename_out '.set'];
    end
    disp(['Loading subject ' filename_in ]);
    EEG_IN = pop_loadset('filename',filename_in,'filepath',filepath);
    EEG_OUT = pop_loadset('filename',filename_out,'filepath',filepath);
    EEG_OUT.icawinv = EEG_IN.icawinv;
    EEG_OUT.icasphere = EEG_IN.icasphere;
    EEG_OUT.icaweights = EEG_IN.icaweights;
    EEG_OUT.icachansind = EEG_IN.icachansind;
    pop_saveset(EEG_OUT, 'filename',['ICA_' filename_out],'filepath',filepath);   
end
