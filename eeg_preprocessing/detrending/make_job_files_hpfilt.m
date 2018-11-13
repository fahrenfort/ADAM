cd ~/lisa_mount/temporient/followup/code/

% cd ~/temporient/followup/code/
lnodes = 1;
cores = 16;
ppn = 15;
mem = 'mem64gb';
lwalltime='2:00:00';

folz = dir('../eegdat/raw/*pp*');
folz([folz.isdir]==0)=[];

scriptname = 'eegpreproc_pop_firws';
cutoff = '0.5';

%-%
subs_per_job = 5;
how_many_jobfiles = ceil(length(folz)/subs_per_job);
how_many_already = length(dir('jobs/*job*.txt'));

k=1;
for filei=1:how_many_jobfiles
    jobid = fopen(['jobs/job' num2str(filei+how_many_already) '.txt'],'w');
    fprintf(jobid, '#PBS -S /bin/bash\n');
    fprintf(jobid, '#PBS -lnodes=%i:cores%i:ppn=%i:%s -lwalltime=%s\n',lnodes,cores,ppn,mem,lwalltime);
    fprintf(jobid,'\nmodule load mcr/v718\nexport MCR_CACHE_ROOT=`mktemp -d /scratch/mcr.XXXXXXXXXX`\n');
    
    if filei<how_many_jobfiles
        q=k+subs_per_job-1;
    else
        q=length(folz);
    end
    
    for ppi=k:q;
        fprintf(jobid,'cp -r $HOME/temporient/followup/eegdat/raw/%s "$TMPDIR" &\n', folz(ppi).name);
    end
    fprintf(jobid,'wait\n');
    for ppi=k:q;
        fprintf(jobid,'$HOME/temporient/followup/code/%s ""$TMPDIR"/%s" "%s" &\n',scriptname,folz(ppi).name,cutoff);
    end
    fprintf(jobid,'wait\n');
    for ppi=k:q;
        fprintf(jobid,'cp "$TMPDIR"/%s/*.set $HOME/temporient/followup/eegdat/processed &\n',folz(ppi).name);
        fprintf(jobid,'cp "$TMPDIR"/%s/*.fdt $HOME/temporient/followup/eegdat/processed &\n',folz(ppi).name);
    end
    
    fprintf(jobid,'wait\nmodule unload matlab/2012b\nmodule unload mcr/v718');
    fclose(jobid);
    k=k+subs_per_job;
end

% --> creates below job textfile format:

% #PBS -S /bin/bash
% #PBS -lnodes=1:cores16:ppn=15:mem64gb -lwalltime=10:00:00
% 
% cp -r $HOME/eegeye/data/pp07 "$TMPDIR"
% cp -r $HOME/eegeye/data/pp08 "$TMPDIR"
% cp -r $HOME/eegeye/data/pp09 "$TMPDIR"
% 
% module load matlab/2012b
% module load mcr/v718
% export MCR_CACHE_ROOT=`mktemp -d /scratch/mcr.XXXXXXXXXX`
% 
% (
%     cd "$TMPDIR"/pp07
%     matlab < run_mcc.m
%     ./tf_conditions_wispc
%     cp *tfdecomp*.mat $HOME/eegeye/data/results
% ) &
% (
%     cd "$TMPDIR"/pp08
%     matlab < run_mcc.m
%     ./tf_conditions_wispc
%     cp *tfdecomp*.mat $HOME/eegeye/data/results
% ) &
% (
%     cd "$TMPDIR"/pp09
%     matlab < run_mcc.m
%     ./tf_conditions_wispc
%     cp *tfdecomp*.mat $HOME/eegeye/data/results
% ) &
% wait
% 
% module unload matlab/2012b
% module unload mcr/v718


%%
fileid = fopen(['bash_' scriptname],'w');
fprintf(fileid,'#!/bin/bash\ncd jobs\n');
for filei=1:how_many_jobfiles
    fprintf(fileid,'qsub job%s.txt\n',num2str(filei+how_many_already));
    fprintf(fileid,'echo qsub job%s.txt\n',num2str(filei+how_many_already));
end
fprintf(fileid,'cd ..');
fclose(fileid);

