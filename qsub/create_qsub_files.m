function create_qsub_files(path_on_lisa, function_name, qsettings, varargin)
% create_qsub_files creates qsub files and bash files to submit qsub jobs for all subjects 
% path_on_lisa contains the path where the function is on lisa e.g. '$HOME/Dropbox/mvpa_scripts'
% function_name contains the function name, e.g. 'classify_RAW_from_eeglab_data'
%
% qsettings contains a bunch of info about how to run the job:
% qsettings.walltime -> how long the job may run
% qsettings.lnodes -> how many nodes
% qsettings.maxcores -> maximum number of cores the script may use (default: 15)
% qsettings.cores -> number of cores the node should have (default: 16)
% qsettings.qsubdir -> where to write the qsub job
%       special value: if qsettings.qsubdir exists, the function uses the path specified in
%       varargin(3), where it replaces $HOME by that the path in qsubdir)
% qsettings.keep_together -> true/false (default: false): determines whether multi-element arguments
% are enumerated together, or whether they are exhaustively combined. This is useful when for
% example different subjects have different condition labels. In that case you can make equally
% sized cell arrays of input files and a cell arrays of condition labels for training/testing.
% qsettings.memory -> memory required (default: none)
% qsettings.repeat -> specify how often each combination of arguments should be repeated, this can
% be useful when running (random permutations). By default, each combination of arguments is
% executed only once.
% qsettings.use_scratch = true (default) -> specify whether to copy files over to the scratch disk
% or not
% qsettings.send_mail = true (default: false) -> the system will forward an e-mail to the e-mail set
% in the $HOME/.forward file when the job starts and finishes
%
% Sometimes Matlab produces errors of compiled versions of the functions during runtime due to the
% fact that some libraries cannot be loaded dynamically. In this case you may get an error like:
% "Can't load '/sara/sw/mcr-r2016a/v901/bin/glnxa64/libmwosgserver.so': dlopen: cannot load any more
% object with static TLS" This can be solved by manually adding the library before runtime. The qsub
% script can take care of this for you by setting: qsettings.preload =
% '/sara/sw/mcr-r2016a/v901/bin/glnxa64/libmwosgserver.so'
% 
% varargin is a series of arguments the function may take each element can be a cell array. All
% combinations of elements are combined, such that each combination becomes a seperate line in your
% qsub file.
%
% By default, if the first first and third arguments of varargin are strings and the second argument
% is a cell array of strings, the script assumes that these contain (1) the input directory on the
% home folder, (2) the input file name(s), and (3) the output directory on the home folder (in that
% order). It will use these to copy the files that the function uses to the scratch space and copy
% the results back to the home folder when done. If the input filenames (2) is specified as a cell
% array of names of files (or combinations of files). If just single files, wildcards * and ? can be
% used, e.g. varargin{2} = '*.set' will try to locate all the .set files and fill up the filenames
% with all .set files found in the input directory.
%
% Example usage:
% create_qsub_files('$HOME/Dropbox/matlab_scripts/eeg_mvpa','classify_RAW_FT_data',qsettings,'/My_experiment/EEGLAB_DATA',{'subj1' 'subj2' 'subj3'},'/My_experiment/MVPA_RESULTS',4,1:4,'linear',true,'1,2,3','4,5,6','7,8,9','10,11,12');
%
% Internal function of the ADAM toolbox by J.J.Fahrenfort, VU, 2014, 2015, 2016, 2018
%
% See also: adam_MVPA_firstlevel

% some default settings
walltime = '23:59:59';
lnodes = '1';
mem = [];
repeat = 1;
use_scratch = true;
keep_together = false;
send_mail = false;
maxcores = [];
cores = 16;
mcr_version = '/v718'; % e.g. '/v718' for matlab 2012b; (or leave empty for the default mcr / latest version of matlab on your machine)
mcr_set_cache = true; % may have to set to true for older versions of MCR, e.g. matlab 2012b
mcr_cache_verbose = false; % produce info about cache creation (can be useful in case of problems)
preload = ''; % e.g. '/sara/sw/mcr-r2016a/v901/bin/glnxa64/libmwosgserver.so'; (if you run into problems loading libraries during runtime)
qsubdir = '';
dotindx = min([find(function_name=='.',1)-1 numel(function_name)]);
bashfilename = regexprep(function_name(1:dotindx),' ','_'); % create a clean file name

% unpack passed settings
v2struct(qsettings);
if ~isempty(mem)
    mem = [':mem' mem];
end
if isempty(cores)
    cores = 12; % maximum nr of jobs to start, always take one less than the number of cores on the node (or even less if you require more memory)
end
if isempty(maxcores)
    maxcores = cores - 1;
end

% locate output path on server through local path
if ~isempty(qsubdir)
    home = strrep(varargin{3},'$HOME',qsubdir);
    if ~exist(home,'dir')
        mkdir(home);
    end
else
    error('you should probably set qsettings.qsubdir (the path on which the remote server is mounted');
end
if ~exist(home,'dir')
    yn = input(['The folder ' home ' does not exist. Do you want to create the file in the present working directory instead (y/n)? '],'s');
    if strcmp(yn,'y')
        home = pwd;
    else
        return
    end     
end

% remove file separator at the end if accidentally present
if strcmp(home(end),filesep)
    home = home(1:end-1);
end

% create job settings
corestxt = [':cores' num2str(cores)];
ppn = [':ppn=' num2str(maxcores)];

% obtain filenames if not already supplied as argument
file = varargin{2};
if ischar(file) && isempty(strfind(file,','))
    if (~isempty(strfind(file,'*')) || ~isempty(strfind(file,'?')))
        if ~isempty(qsubdir)
            file = [strrep(varargin{1},'$HOME',qsubdir) filesep file];
        else
            file = [varargin{1} filesep file];
        end
        file = dir(file);
        files = {file(:).name};
        if isempty(files)
            error('cannot find any files in that location and/or that match that wildcard pattern');
        else
            varargin{2} = files;
        end
    end
end

% specify input/output on scratch space
timeout = '';
indir =  varargin{1};
outdir = varargin{3};
if use_scratch
    timeout = 'timeout $timetorun ';
    scratchindir = '"$TMPDIR"/input';
    scratchoutdir = '"$TMPDIR"/output';
    scratchlogdir = '"$TMPDIR"/matlablogdir';
    varargin{1} = scratchindir;
    varargin{3} = scratchoutdir;
end

% obtain argument combinations to create a loop around each var combination in varargin
if keep_together
    combMat = enumcombs(varargin{:});
else
    combMat = allcombs(varargin{:});
end

% multiply the set by repeat (e.g. when doing many random permutations)
combMat = repmat(combMat,repeat,1);

% when doing iterations, create a separate qsub file for every subject
if repeat > 1
    allfiles = combMat(:,2);
    uniquefiles = unique(allfiles);
    allMat = cell(size(uniquefiles));
    for cFiles = 1:numel(uniquefiles)
        allMat{cFiles} = combMat(strcmp(uniquefiles{cFiles},allfiles),:);
    end
else % or run as normal
    allMat = {combMat};
end

% create qsub job for each subject (in case repeat > 1) and a bash file for all jobs 
qsubfiles = {};
for cMat = 1:numel(allMat)
    combMat = allMat{cMat};
    for cQsubs = 1:size(combMat,1)
        if (mod(cQsubs,maxcores) == 1 && repeat == 1) || maxcores == 1 || cQsubs == 1
            % initialize
            copyin = [];
            copyout = [];
            line = [];
            % create qsub file to submit all subjects, add trailing nr if it already exists
            c = 1;
            qsubfile = sprintf([home filesep 'qsub_' bashfilename '_%03d'], c);
            while exist(qsubfile,'file')
                c = c + 1;
                qsubfile = sprintf([home filesep 'qsub_' bashfilename '_%03d'], c);
            end
            qsubfiles{end+1} = qsubfile;
            fout = fopen(qsubfile,'w');
            fprintf(fout,'#PBS -S /bin/bash\n');
            fprintf(fout,['#PBS -lnodes=' lnodes corestxt ppn mem ' -lwalltime=' walltime '\n']);
            fprintf(fout,'echo "Job $PBS_JOBID started at `date`"');
            if send_mail
                fprintf(fout,' | mail $USER -s "Job $PBS_JOBID"');
            end
            fprintf(fout,'\n');
            if use_scratch
                % build in some time to copy output back by breaking off jobs before the end
                fprintf(fout,'module load sara-batch-resources\n');
                % make sure directories exists
                fprintf(fout,['mkdir -p ' scratchindir ' &\n']);
                fprintf(fout,['mkdir -p ' scratchoutdir ' &\n']);
                fprintf(fout,['mkdir -p ' scratchlogdir ' &\n']);
                % build in some time to copy output back by breaking off jobs before the end
                fprintf(fout,'(( timetorun = $SARA_BATCH_WALLTIME - 600 ))\n');
            end
            % load matlab runtime
            fprintf(fout,['module load mcr' mcr_version '\n']);
            if mcr_set_cache
                fprintf(fout,'export MCR_CACHE_ROOT=`mktemp -d "$TMPDIR"/mcr.XXXXXXXXXX`\n');
            end
            if mcr_cache_verbose
                fprintf(fout,'export MCR_CACHE_VERBOSE=true\n');
            end
            if ~isempty(preload)
                fprintf(fout,['export LD_PRELOAD=' preload  '\n']);
            end
            if use_scratch
                fprintf(fout,['export MATLAB_LOG_DIR=' scratchlogdir '\n']);
            end
            % make sure directories exists
            fprintf(fout,['mkdir -p ' outdir ' &\nwait\n']);
        end
        % copy all the files to and from scratch
        if use_scratch
            filename = combMat{cQsubs,2};
            filenames = regexp(filename,',','split');
            for c = 1:numel(filenames)
                filename = filenames{c};
                if regexp(filename(end-3:end),'\....') == 1
                    filename = filename(1:end-4);
                end
                if isempty(strfind(copyin,['cp -u ' indir '/' filename '.* ' scratchindir]))
                    copyin = [copyin 'cp -u ' indir '/' filename '.* ' scratchindir ' &\n'];
                end
            end
        end
        % commands to issue in qsub file
        line = [ line timeout path_on_lisa '/' function_name];
        for cArgs = 1:size(combMat,2)
            line = [line ' "' combMat{cQsubs,cArgs} '" '];
        end
        line = [line ' &\n'];
        % close qsub file once all cores are used or all commands have been issued
        if (mod(cQsubs,maxcores) == 0 && repeat > 1 && ~(cQsubs==size(combMat,1)))
            % pause till all previous are done
            line = [line 'wait\n'];
        elseif (mod(cQsubs,maxcores) == 0 && repeat == 1) || maxcores == 1 || cQsubs == size(combMat,1)
            % copy to scratch
            if use_scratch
                copyin = [copyin 'wait\n'];
                fprintf(fout,copyin);
            end
            % write individual commands to qsub job
            line = [line 'wait\n'];
            fprintf(fout,line);
            % copy back to home
            if use_scratch
                copyout = ['cp -u -r ' scratchoutdir '/* ' outdir ' &\nwait\n'];
                fprintf(fout,copyout);
            end
            if send_mail
                fprintf(fout,'echo "Job $PBS_JOBID finished at `date`" | mail $USER -s "Job $PBS_JOBID"\n');
            else
                fprintf(fout,'echo "Job $PBS_JOBID finished at `date`"\n');
            end
            fclose(fout);
            fclose('all');
        end
    end
end

% create a bash file to qsub all the generated qsubs
c = 1;
bashfile = sprintf([home filesep '..' filesep 'bash_' bashfilename '_%03d'], c);
while exist(bashfile,'file')
    c = c + 1;
    bashfile = sprintf([home filesep '..' filesep 'bash_' bashfilename '_%03d'], c);
end
fout = fopen(bashfile,'w');
fprintf(fout,'#!/bin/bash\n');
fprintf(fout,['cd ' home(max(strfind(home,'/'))+1:end) '\n']);
for cqsubs = 1:numel(qsubfiles)
    do = ['qsub ' qsubfiles{cqsubs}(max(strfind(qsubfiles{cqsubs},'/'))+1:end) '\n'];
    fprintf(fout,do);
    fprintf(fout,['echo ' do]);
end
fprintf(fout,'cd ..\n');
fclose(fout);
fclose('all');
disp(['To submit all jobs bash the following file: ' bashfile]);

return

function combMat = allcombs(varargin)
% make all possible combinations of all arguments
for c = 1:nargin
    if isnumeric(varargin{c})
        varargin{c} = strsplit(num2str(varargin{c}),' ');
    end
    if ischar(varargin{c})
        varargin{c} = {varargin{c}};
    end
end
sizeVec = cellfun('prodofsize', varargin);
indices = fliplr(arrayfun(@(n) {1:n}, sizeVec));
[indices{:}] = ndgrid(indices{:});
combMat = cellfun(@(c,i) {reshape(c(i(:)), [], 1)}, ...
    varargin, fliplr(indices));
combMat = [combMat{:}];
return

function enumMat = enumcombs(varargin)
% enumerate arguments, keeping the elements of each argument together with
% other elements
for c = 1:nargin
    if isnumeric(varargin{c})
        varargin{c} = strsplit(num2str(varargin{c}),' ');
    end
    if ischar(varargin{c})
        varargin{c} = {varargin{c}};
    end
end
maxArg = 1;
for c = 1:nargin
    if numel(varargin{c}) > maxArg;
        maxArg = numel(varargin{c});
    elseif numel(varargin{c}) ~= maxArg && numel(varargin{c}) ~= 1
        disp(numel(varargin{c}));
        error('all argument enumerations larger than 1 should be the same size, or I cannot keep arguments together. It might be that you are trying to run with qsettings.keep_together = true in combination with more than one channelset in one go.');
    end
end
for cQsub = 1:maxArg
    for cArg = 1:nargin
        if numel(varargin{cArg}) > 1
            enumMat(cQsub,cArg) =  varargin{cArg}(cQsub);
        else
            enumMat(cQsub,cArg) = varargin{cArg};
        end
    end
end
return