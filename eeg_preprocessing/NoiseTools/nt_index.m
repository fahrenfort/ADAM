function [status,p]=nt_index(name,p,forceUpdate)
%[status,p]=nt_index(name,p,forceUpdate) - index data files & directories
%
%  status: 1: needed indexing, 0: didn't, -1: failed
%  p: parameter structure
%
%  name: name(s) of file(s) or directory to index
%  p: parameters
%  forceUpdate: if true force indexing [default: false]
%
% NoiseTools
nt_greetings;

if nargin<3 || isempty(forceUpdate); forceUpdate=0; end 
if nargin<2||isempty(p) % set default parameters
    p=[];
    p.scale=1000;
    if nargin >= 1; p.name=name; end    
end
if nargin<1 || isempty(name)
    p.name=[];
    status=-1;
    return;  % just return default parameters
end

status=-1; % failed by default
updateFlag=0;   % don't update unless necessary
if forceUpdate; updateFlag=1; end

% parse 'name' into path, etc., check it for various issues
if ~ischar(name); error('name should be a string'); end
avoid=['[',1:31, 127,']'];
if regexp(name,avoid) 
    disp('bad character in file name, skip:'); disp(['   >',name,'<']); 
    return; 
end
if name=='.'; name=pwd; end
if name(end)=='/'; name=name(1:end-1); end 
[PATHSTR,NAME,EXT]=fileparts(name);
if strcmp(EXT,'idx')
    disp(['warning: ', name, ' might be index file']); 
end
if isempty(PATHSTR)             % interpret relative to current directory
    name=[pwd,filesep,name];    % need full path to safely use 'exist'
end
[PATHSTR,NAME,EXT]=fileparts(name); 
if 2==exist(name) 
    d=dir(name);
    filename=d.name;            % --> same case as file system
    PATHSTR=cd(cd(PATHSTR));    % --> same case as file system
    name=[PATHSTR,filesep,filename];
elseif 7==exist(name)
    name=cd(cd(name));          % --> same case as file system
    [PATHSTR,NAME,EXT]=fileparts(name); 
else
    disp(name);
    error('...is neither file nor directory');
end


hhh=[]; % this structure will contain info about this file or directory
iii=[]; % this structure will contain the data index
hhh.name=name;
hhh.time_indexed=now;
hhh.failed=0; % OK by default

% test whether we're processing a file or a directory
if 2==exist(name) 
    hhh.isdir=0;
elseif 7==exist(name)
    hhh.isdir=1;
else
    disp(name);
    error('...is neither file nor directory');
end

% special case:
% CTF data are stored as directory, pretend it's a file
if numel(name>=3) && strcmp(name(end-2:end), '.ds')
    hhh.isdir=0;
end
        
% create an index directory if it doesn't exist
idxDir=[PATHSTR,filesep,'nt_idx'];
if 7 ~= exist(idxDir) 
    disp(['creating index directory ', idxDir]);
    mkdir (idxDir);
    updateFlag=1;
end

% check if there is already an up-to-date index file for 'name'
idxName=[idxDir,filesep,NAME,EXT,'.idx'];
hhh.idxName=idxName;
if ~2==exist(idxName); updateFlag=1; end
if exist(idxName) && (dateModified(idxName) < dateModified(name)) % out of dat
    updateFlag=1;
end

%{
Processing depends on whether 'name' is a file or a directory.
If a file, we calculate statistics to index the data within that file.
If a directory, we aggregate statistics on its files and subdirectories. 
%}

if hhh.isdir % directory 
        
    disp([name,filesep]);
    
    % check that 'name' matches a name listed in parent directory (catch upper/lowercase inconsistencies)
    d=dir(PATHSTR);
    OKflag=0;
    for iFile=1:numel(d)
        if strcmp(d(iFile).name,[NAME,EXT])
            OKflag=1;
        end
    end
    if ~OKflag; error(['''', NAME,EXT, ''' does not match real file name']); end   

    % list items in this directory
    d=dir(name);
    iGood=ones(numel(d),1);
    nskip=0;
    for k=1:numel(iGood)       % weed out irrelevant/bad files
        if strcmp(d(k).name,'.') || strcmp(d(k).name,'..')  % me & parent dirs
            iGood(k)=0; 
        elseif d(k).name(1)=='.'                       	% files starting with '.'
            iGood(k)=0; nskip = nskip+1; 
            disp(['skip, starts with ''.'': ',name,filesep,d(k).name]);
        end                             
        if strcmp(d(k).name,'nt_idx')                 	% index directory
            iGood(k)=0; nskip=nskip+1; 
        end      
        if any(d(k).name<33) || any(d(k).name==127)  	% files with bad chars in names
            iGood(k)=0; nskip=nskip+1;
            disp(['skip, bad char in name: ',name,filesep,d(k).name]);
        end 
        if isempty(d(k).date)                            % no date field ==> invalid file (soft link?)
            iGood(k)=0; nskip=nskip+1; 
            disp(['skip, invalid file: ',name,filesep,d(k).name]);
        end 
    end
    d=d(iGood~=0);
    nfiles=numel(d);
    
    % recursively index files in this directory
    for iFile=1:nfiles
        if d(iFile).name(end)=='/'
            disp(d)
        end
        if 1==nt_index([name,filesep,d(iFile).name],p,forceUpdate) % recurse
            updateFlag=1;  % one of my files updated, update me too
        end        
    end   
    
    % purge my index directory of any orphan index files (no associated data file)
    dd=dir(idxDir);
    iGood2=ones(numel(dd),1);
    for k=1:numel(iGood2)      
        if dd(k).name(1)=='.'; iGood2(k)=0; end             % files starting with '.'
        if strcmp(dd(k).name,'nt_idx'); iGood2(k)=0; end    % index directory
    end
    dd=dd(iGood2~=0);
    iGood3=ones(numel(dd),1);
    for iFile=1:numel(iGood3)
        [~,NAME2,EXT2]=fileparts(dd(iFile).name);           % name of index file
        theFile=[PATHSTR,filesep,NAME2];                    % associate data file
        if 2~=exist(theFile) ...        % neither file...
                && 7~=exist(theFile)    % nor directory...
            disp(['>',dd(iFile).name,'<'])
            disp([theFile, ' not found, ']);
            disp(['deleting orphan index file ',[idxDir,filesep,dd(iFile).name]]);
            delete([idxDir,filesep,dd(iFile).name]);
            iGood3(iFile)=0;
        end
    end
    dd=dd(iGood3~=0);
    
    % all the files in this directory & subdirectories are now checked/updated
    
    % merge info about this directory, files, & subdirectories into index
    if updateFlag
        
        % info about this directory
        hhh.dir=d;                  % my directory structure, excluding bad files
                
        % init aggregated statistics
        hhh.nfiles=uint64(1);   % number of files
        hhh.ndata=0;            % number of data files
        hhh.nbad=0;             % number of bad files
        hhh.nskip=nskip;        % number of files skipped
        hhh.ndirs=uint64(1);    % number of directories
        hhh.bytes=uint64(0);    % number of bytes
        hhh.ntypes=[];          % list with number of files of each type
        hhh.depth=1;            % depth of hierarchy (1=leaf)
        
        % init info variables for each file/directory in this directory
        hhh.filelist.bytes=zeros(nfiles,1,'uint64');       % bytes (file) or total bytes (directory)
        hhh.filelist.isdir=nan(nfiles,1);                  % directory?
        hhh.filelist.nfiles=zeros(nfiles,1,'uint64');      % number of files (including files in subdirectories)
        
        % visit each file or directory, aggregate information
        for iFile=1:nfiles
            
            % get info from this file's index file to save time
            load('-mat',[name,filesep,'nt_idx',filesep,d(iFile).name,'.idx'], 'hh');  

            % info specific to this file/directory
            hhh.filelist.bytes(iFile)=hh.bytes;
            hhh.filelist.isdir(iFile)=hh.isdir;
            hhh.filelist.nfiles(iFile)=hh.nfiles;
            
            % aggregate info
            hhh.bytes=hhh.bytes+hh.bytes;
            hhh.ndirs=hhh.ndirs+hh.ndirs;
            hhh.ndata=hhh.ndata+hh.ndata;
            hhh.nbad=hhh.nbad+hh.nbad;
            hhh.nfiles=hhh.nfiles+hh.nfiles;   
            hhh.bytes=hhh.bytes+hh.bytes;
            hhh.nskip=hhh.nskip+hh.nskip;
            hhh.depth=max(hhh.depth,1+hh.depth);
            
            % aggregate counts of each file type
            types=myfieldnamesr(hh.ntypes);
            for iType=1:numel(types)
                if isfield(hhh.ntypes,types(iType))
                    %eval(['hhh.ntypes.',types{iType},'=hhh.ntypes.',types{iType},'+hh.ntypes.',types{iType},';']);
                    hhh.ntypes.(types{iType})=hhh.ntypes.(types{iType}) + hh.ntypes.(types{iType});
                else
                    eval(['hhh.ntypes.',types{iType},'=hh.ntypes.',types{iType},';']);
%                   hhh.ntypes.(types{iType})= hh.ntypes.(types{iType}); dumb matlab can't do this properly
                end
            end
            
        end
        
        % merge the indexes of files & subdirectories to create an aggregate index  
        
        iii=merge_file_indexes(d,[PATHSTR,filesep,NAME]);
        
    end % if updateflag
        
else  % file
    
    %disp(name)
    hhh.isdata=0;
    
    if numel(name>=3) && strcmp(name(end-2:end), '.ds') % intercept CTF data
        [a,b,c] = fileparts(name);
        name=[name,filesep,b,'.meg4'];
    end
       
    if updateFlag
        
        % info common to all files
        hhh.nfiles=uint64(1); % just me
        d=dir(name);
        hhh.bytes=uint64(d.bytes);
        hhh.sr=[];
        hhh.depth=0;
        
        % default values:
        hhh.ndirs=uint64(0); 
        hhh.nbad=0; 
        hhh.ndata=0;
        hhh.nskip=0;
        
        % determine file type
        [isdata,type]=filetype(name);
        hhh.isdata=isdata;
        hhh.type=type;
        
        % set field in hhh.ntypes
        fixedtype=strrep(type,':','___'); % biosig uses ':' in type names
        try
            eval(['hhh.ntypes.',fixedtype,'=1;']);
        catch
            disp(['hhh.ntypes.',fixedtype,'=1;']);
            disp(name);
            disp(type);
            warning('eval failed');
        end
        
        % data: read it
        if hhh.isdata
            x=[];
            hhh.size=[];   
            hhh.originalsize=[]; % before reshape/transpose
            hhh.ndata=1;
            [a,b,c]=fileparts(type);
            if strcmp(b,'matlab')
                % read matlab variable from file
                variable_name=c(2:end); % was coded as extension to type
                x=readmatlab(name,variable_name);
            elseif strcmp(type,'unknown') || strcmp(type,'matlab_non_numeric')
                % shouldn't happen
                error('!');
            else
                % some data file, try to read with biosig
                try
                    h=sopen(name);
                    hhh.sr=h.SampleRate;
                catch ME
                    hhh.failed=1;
                    disp(name);
                    warning('...sopen failed');
                    disp(ME);
                end
                try
                    x=sread(h);
                catch ME
                    hhh.failed=1;
                    disp(name)
                    disp(ME);
                    warning('...sread failed');
                    x=sread(h);
                end
                sclose(h);
            end
            
            % transpose if appropriate (this is a kludge)
            hhh.originalsize=size(x);
            if ndims(x)>2
                % more than 3 dims, merge last dimensions
                sizes=size(x);
                x=reshape(x,prod(sizes(1:end-1)),sizes(end));
                disp(['reshape -->', num2str(size(x))]);
            end
            if size(x,1)<size(x,2) 
                % wider than tall: transpose
                x=x'; 
                disp(['transpose --> ',num2str(size(x))]);
            end
            hhh.size=size(x);
            nt_whoss;

            if ~isempty(x)
                % calculate index
                dsratio=100;
                iii.card=[]; iii.min=[]; iii.max=[]; iii.mean=[]; iii.ssq=[];
                iii=nt_idx(x,dsratio,iii);
            end % else iii==[]
        end
    end
end   

if updateFlag
    status=1;
    hh=hhh; ii=iii;
    save(idxName, 'hh','ii');
    disp(idxName)
else 
    status=0;
end
end % function [status,p]=nt_index(name,p,forceUpdate)

function ii=index(x,p)
% index data
if ndims(x)>2; error('!'); end
[ii.nsamples,ii.nchans]=size(x);
ii.scale=p.scale;
ii.p=p;
npairs=floor(ii.nsamples/p.scale);
size(x)
x_extra=x(npairs*p.scale+1:end,:);
x=x(1:npairs*p.scale,:);
x=reshape(x,[p.scale,npairs,ii.nchans]);
ii.min=squeeze(min(x,[],1))';
ii.max=squeeze(max(x,[],1))';
if ~isempty(x_extra)
    [size(ii.min) size(x_extra)]
    ii.min=[ii.min;min(x_extra,[],1)];
    ii.max=[ii.max;max(x_extra,[],1)];
end
end

function date=dateModified(name)
% modification date of file or directory
[PATHSTR,NAME,EXT]=fileparts(name);
if isempty(PATHSTR); error('!'); end
date=[];
if 2==exist(name) % I'm a file, I own my date.
    d=dir(name); % get directly from file
    date=d.datenum;
elseif 7==exist(name) % I'm a directory, my parent own's my date
    d=dir(PATHSTR);
    for iFile=1:numel(d)
        %disp(d(iFile).name)
        if strcmp(d(iFile).name,[NAME,EXT])
            date=d(iFile).datenum; % get indirectly from parent directory
            break
        end
    end
else
    disp(name)
    error('!');
end
if isempty(date) 
    disp(['>',name,'<']);
    error('!'); 
end
end % function date=dateModified(name)

        

function [isdata,type]=filetype(name)
% try to guess type and whether it's data
EXTENSIONS_TO_SKIP={'.idx', '.zip','.txt','.pdf','.doc','.docx','.ppt','.pptx','.xls','.html','.rtf',...
    '.jpg', '.png', '.tif','.tiff','.js', '.md', '.m', '.py', '.rar', '.wav', '.eps', '.pdfsync',...
    '.avi', '.PDF', '.gz', '.zip'};
[PATHSTR,NAME,EXT]=fileparts(name);
isdata=0; type='unknown'; transpose=0; % default
d=dir(name);
if d.bytes==0 
    type='empty';
elseif ~isempty(EXT) && any(strcmpi(EXT,EXTENSIONS_TO_SKIP))
    isdata=0; type=lower(EXT); type=type(2:end); % intercept common types
    disp(['skip (extension): ',name])
else
    fid=fopen(name);
    firstbytes=fread(fid,8,'uchar');
    fclose(fid);
    if ~isempty(EXT) && strcmp(EXT,'.mat') || (numel(firstbytes)>=4 && all(firstbytes(1:4)'=='MATL')) 
        % matlab file
        try
            s=whos('-file',name);
        catch ME
            disp('name');
            disp('... whos failed');
            disp(ME)
            type=[]; return
        end
        % find which variables are numeric
        numerics={'double','single','int64','int32','int16','int8'};
        matrix=strcmp(repmat({s.class},numel(numerics),1), ...
            repmat(numerics',1,numel(s))); 
%         idx=find(any(matrix));
%         if isempty(idx)
        if ~any(matrix)
            % no numeric variables
            isdata=0; type='matlab_non_numeric';
        else
            % some variables are numeric, choose the biggest one
            sizes=zeros(numel(s),1);
            for iVariable=1:numel(s)
                if any(strcmp(s(iVariable),numerics))
                    sizes(iVariable)=prod(s(iVariable).size);
                else 
                    sizes(iVariable)=0;
                end
            end
            [~,biggest]=max(prod(sizes));
            isdata=1; type=['matlab.',s(biggest).name];
            disp(name);
            disp(['mat file, multiple numeric variables, chosing: ''', s(biggest).name, ''', size:',num2str(s(biggest).size)]);
        end
%     elseif strcmp(char(firstbytes(2:8))','BIOSEMI')
%         isdata=1; type='biosemi_bdf';
    else    % hand over to biosig
        try
            h=sopen(name);
            type=h.TYPE;
            sclose(h);
        catch ME
            disp(name);
            warning('... sopen failed');
            disp(ME);
            type='unknown'; return
        end                
        if strcmp(type,'unknown')
            isdata=0; 
        else 
            isdata=1;
        end
    end
end
end % function [isdata,type]=filetype(name)
    
function x=readmatlab(name,varname)
% read varname from matlab file
load('-mat',name,varname);
eval(['x=',varname, ';']);
end % function x=readmatlab(name,varname)


% recursive tally of field names at all depths
function s=myfieldnamesr(x)
if ~isstruct(x); s=[]; return; end
fields=fieldnames(x);
s={};
for iField=1:numel(fields)
    xx=getfield(x,fields{iField});
    if isa(xx,'struct')
        subfields=myfieldnamesr(xx);
        for iSubfield=1:numel(subfields)
            s=[s,[char(fields(iField)),'.',char(subfields(iSubfield))]];
        end
    else
        s=[s,fields{iField}];
    end
end
end % function x=readmatlab(name,varname)


% visit index files of all files in 'd', merge into aggregate index
function iii=merge_file_indexes(d, dname)
    iii=[];
    % load all indexes into cell array
    all_indexes=[];       
    nfiles=numel(d);
    for iFile=1:nfiles
        name=d(iFile).name;
        load('-mat',[dname,filesep,'nt_idx',filesep,d(iFile).name,'.idx'], 'ii'); 
        all_indexes{iFile}=ii;
    end
    % estimate total size of statistics
    statNrows=[]; % cell array of stat sizes
    statNcols=[]; % cell array of stat sizes
    for iFile=1:nfiles
        ii=all_indexes{iFile};
        if isempty(ii); continue; end
        statNames=fieldnames(ii);
        for iField=1:numel(statNames)
            if ~isfield(statNcols,statNames{iField})
                %setfield(statSizes,statNames{iField})=0;
                statNrows.(statNames{iField})=0;
                statNcols.(statNames{iField})=0;
            end
            [nrows,ncols]=size(getfield(ii, statNames{iField}));
            %tmp=getfield(statSizes,statNames{iField});
            %setfield(statSizes,statNames{iField},tmp);
            statNrows.(statNames{iField})=statNrows.(statNames{iField})+nrows;
            statNcols.(statNames{iField})=max(statNcols.(statNames{iField}),ncols);
        end
    end
    if isempty(statNcols)
        return; % no indexes to summarize
    end
    statNames=fieldnames(statNcols);
    iCounter=[];
    for iName=1:numel(statNames) 
        statName=statNames{iName};
        %setfield(iii,statNames{iName},zeros(statSizes(iName)));
        iii.(statNames{iName})=zeros(statNrows.(statName),statNcols.(statName));
        %setfield(iCounter, statnames{iName}, 0); 
        iCounter.(statNames{iName})=0;
    end 
    for iFile=1:nfiles
        ii=all_indexes{iFile};
        for iName=1:numel(statNames)
            statName=statNames{iName};
            if isfield(ii,statName)
                tmp=ii.(statName);
                offset=iCounter.(statName);
                iii.(statName)(offset+(1:size(tmp,1)),1:size(tmp,2))=tmp;
                iCounter.(statName)=iCounter.(statName)+size(tmp,1);
            end
        end
    end
end % function iii=merge_file_indexes(d)        

           
        
        
            
   