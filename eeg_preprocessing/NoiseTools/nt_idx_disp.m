function nt_idx_disp(name,field,explainflag)
%nt_idx_disp(name,field,explainflag) - display contents of index file
%
%  name: name of file
%  field: field to display in detail
%  explainflag: if true, explain contents of each field
%
% NoiseTools
nt_greetings;

if nargin<1; name=pwd; end
if nargin<2; field=[]; end
if nargin<3; explainflag=0; end

% check name, parse into path, etc.
if ~ischar(name); error('name should be a string'); end
avoid=['[',1:31, 127,']'];
if regexp(name,avoid); 
    disp('bad character in file name, skip:'); disp(['   >',name,'<']); 
    return; 
end
if name=='.'; name=pwd; end
if name(end)=='/'; name=name(1:end-1); end % remove trailing slash
[PATHSTR,NAME,EXT]=fileparts(name);
if strcmp(EXT,'idx'); 
    disp(['warning: ', name, ' might be index file']); 
end
if isempty(PATHSTR); % interpret relative to current directory
    name=[pwd,filesep,name]; % full path, safe to use 'exist'
end
[PATHSTR,NAME,EXT]=fileparts(name); 
if 2==exist(name) 
    d=dir(name);
    filename=d.name;            % match case to file system
    PATHSTR=cd(cd(PATHSTR));    % match case to file system
    name=[PATHSTR,filesep,filename];
elseif 7==exist(name)
    name=cd(cd(name));          % match case to file system
    [PATHSTR,NAME,EXT]=fileparts(name); 
else
    disp(name);
    error('...is neither file nor directory');
end

idxDir=[PATHSTR,filesep,'nt_idx'];
idxFile=[idxDir,filesep,NAME,EXT,'.idx'];

if 2~=exist(idxFile); 
    disp(idxFile)
    disp('index file not found');
    return
end

load('-mat',idxFile);  % loads hh, ii
a.hh=hh;
a.ii=ii;

if ~isempty(field)
    eval(['a=a.',field,';']);
end

if ~explainflag
    if isnumeric(a)
        figure(100); clf; plot(nt_demean(a)); 
    else
        disp(a);
    end
elseif ~isempty(a)
    fieldNames=fieldnames(a);
    for iFieldName=1:numel(fieldNames)
        fieldName=fieldNames{iFieldName};
        switch fieldName
            case 'hh'
                disp('hh: header and metadata');
            case 'ii'
                disp('ii: index describing the data')
            case 'name'
                disp('name: name of file or directory that was indexed');
            case 'idxName'
                disp('idxName: name of index file');
            case 'time_indexed'
                disp('time_indexed: date/time at which this was indexed (datenumber)');
            case 'failed'
                disp('failed: 1 = indexing failed');
            case 'isdir'
                disp('isdir: 1 = directory');
            case 'dir'
                disp('dir: directory list (for directory) or entry (for file)') 
            case 'filelist'
                disp('filelist: info for each file in this directory');
            case 'isdata'
                disp('isdata: 1 = file was recognized as containing indexable data');
            case 'nfiles'
                disp('nfiles: total number of indexed files within this directory & subdirectories');
            case 'bytes'
                disp('bytes: total number of bytes within this file or directory & subdirectories');
            case 'date'
                disp('date: directory date field for this file or directory');
            case 'sr'
                disp('sr: sampling rate, if known');
            case 'depth'
                disp('depth: depth of the subdirectory hierarchy');
            case 'ndirs'
                disp('ndirs: total number of directories within this directory & subdirectories');
            case 'nbad'
                disp('nbad: dunno what this means...');
            case 'ndata'
                disp('ndata: total number of data files in this directory & subdirectories');
            case 'nskip'
                disp('nskip: total number of files that were skipped in this directory & subdirectories');
            case 'ext'
                disp('ext: extension of this file''s name');
            case 'type'
                disp('type: type of this data file');
            case 'ntypes'
                disp('ntypes: number of files of each type in this directory & subdirectories');
            case 'size'
                disp('size: dimensions of matrix being indexed')
            case 'originalsize'
                disp('originalsize: original dimensions before transposing');
            case 'min'
                disp('min: minimum over interval of samples');
            case 'max'
                disp('max: maximum over interval of samples');
            case 'mean'
                disp('mean: mean over interval of samples');
            case 'var'
                disp('var: variance over interval of samples')
            case 'card'
                disp('card: cardinality of interval of samples');
            otherwise
                disp([fieldName,': ?']);
        end
    end
end
              
                
