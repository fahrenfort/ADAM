function [p,data]=nt_read_data(fname,flag)
%[p,data]=nt_read_data(fname,flag) - read data from file
%
% 
%  fname: file to read
%  flag: specify how to deal with complex data
%       0: walk through and select 
%       1: choose first of each [default]
%       2: return all as struct
%  
% 
VERBOSE=1;

if nargin < 1 ; error('!'); end
if nargin < 2; flag=1; end

if ~isa(fname, 'char'); 
    error('filename is not char string');
end
if exist(fname,'file')~=2 
    error(['file >', fname, '< not found']);
end

% standard fields
p.fname=fname;
p.read_with=[];
p.sr=[];

% intercept directories
if exist(fname,'dir') 
    dname=fname; clear fname;
    if VERBOSE; disp('directory'); end    
    d=dir(dname);
    fnames=char(d.name);
    idx=find(fnames(:,1)~='.');  % remove '.' and '..' and invisible files
    d=d(idx);
    fnames=fnames(idx,:);
    if numel(d)==0
        error(['directory >',fname,'< is empty']);
    end
    
    
    % separate directories and files
    didx=find([d.isdir]);
    fidx=find(~[d.isdir]);
    fnames=fnames([didx, fidx],:);
    
    switch flag
        case 2
            % return all as array
            data={}; p={};
            for iFile=1:size(fnames,1);
                [pp,dd]=nt_read_data([dname,deblank(fnames(iFile,:))],flag);
                p{iFile}=pp;
                data{iFile}=dd;
            end
            return
        case 1
            % return first
            if VERBOSE; disp(['choosing file: ', fnames(1,:)]); end
            [p,data]=nt_read_data([dname,deblank(fnames(1,:))],flag);
            return
        case 0
            % choose one
            ;
    end
    
   % count files within the directories
    nfiles=zeros(numel(didx),1);
    for k=1:numel(didx)
        dd=dir([data,'/',d(didx(k)).name]);
        fns=char(dd.name);
        idx=find(fns(:,1)~='.');  % remove '.' and '..' and invisible files
        nfiles(k)=numel(idx);
    end
    
    % size of the files
    mbytes=[d(fidx).bytes]'/1024;
   
    % string arrays to put in dialog list
    a=repmat(' (', numel(d),1);
    if numel(didx)>0
        b=cellstr(num2str(nfiles, '%9d'));
    else
        b=[]; % stupid matlab!
    end
    if numel(fidx)>0
        b=[b;cellstr(num2str(mbytes,'%0.1f'))];
    end
    b=char(b);
    c=[repmat(' files)', numel(didx),1); repmat(' Mb)   ', numel(fidx),1)];
     
    % which directory or file is user interested in?    
    i=listdlg('liststring',cellstr([fnames,a,b,c]),...
        'name', 'Select file:', ...
        'listsize', [300 300], ...
        'OKstring','Select',...
        'PromptString','choose file');
    
    if numel(i)==1; 
        [p,data]=nt_read_data([dname,deblank(fnames(i,:))]); 
    elseif isempty(i)
        p=[]; data=[];
    else
        % load a bunch of files
        p.fnames=fnames(i);
        for iFile=1:numel(i);
            [pp,dd]=nt_read_data([dname,deblank(fnames(i(iFile),:))]);
            p.p{iFile}=pp;
            data{iFile}=dd;
        end
    end
    return
 
end

% intercept .mat files
if numel(fname)>4 & fname(end-3:end)=='.mat'
    p.read_with='.mat file';
    if VERBOSE; disp('mat file'); end
    % list variables in file, ask user to choose one
    switch flag
        case 2
            % load all as struct
            data=load(fname);
            return
        case 1
            % load first variable in mat file
            S=whos('-file',fname);
            if VERBOSE; disp(['choosing variable: ',S(1).name]); end
            data=load(fname,deblank(S(1).name));
            while isstruct(data);
                % load first field in structure
                if numel(data)>1
                    data=data(1); 
                end
                S=fieldnames(data);
                if VERBOSE; disp (['choosing field: ',S{1}]); end
                data=getfield(data,S{1});
            end
            return
        case 0
            % choose
    end
    S=whos('-file',fname);
    var_names=char(S.name);
    var_sizes=round([S.bytes]/1024)';
    a=repmat(' (', size(var_names,1),1);
    b=cellstr(num2str(var_sizes, '%9d'));
    b=char(b);
    c=[repmat(' Mb)', size(var_names,1),1)];
    i=listdlg('liststring',cellstr([var_names,a,b,c]),...
        'name', ['Select variable in file ',fname], ...
        'listsize', [600 300], ...
        'OKstring','Select',...
        'PromptString','select:');
    if isempty(i); data=[]; return; end
    if nargout>1;
        data=load(fname,deblank(var_names(i,:)));
        % if it's a structure, list fields, ask user to choose one
        while isstruct(data);
            if numel(data)>1
                i=listdlg('liststring',cellstr(S),...
                    'name', ['Select element of strucure array ',var_names(i,:)], ...
                    'listsize', [600 300], ...
                    'OKstring','Select',...
                    'PromptString','select:');
                if i ; data=data(i); end
            end
            S=fieldnames(data);
            i=listdlg('liststring',cellstr(S),...
                'name', ['Select field in struct ',var_names(i,:)], ...
                'listsize', [600 300], ...
                'OKstring','Select',...
                'PromptString','select:');
            if i ; data=getfield(data,S{i}); end
        end
    end
    return
end

% intercept Yokogawa files
if numel(fname)>4 & (fname(end-3:end)=='.con' | fname(end-3:end)=='.sqd')
    p.read_with='yokogawa 2013';
    p.acq_cond = getYkgwHdrAcqCond(fname);
    p.channel_info=getYkgwHdrChannel(fname);
    p.system_info=getYkgwHdrSystem(fname);
    p.event=getYkgwHdrEvent(fname);
    % read other info?
    p.sr=p.acq_cond.sample_rate;
    if nargout>1;
        data=getYkgwData(fname)';
    end
    return
end
   
         
% select file reader among those available
has_ft_reader=0; 
has_sopen=0;
if 2==exist('ft_read_header');
    has_ft_read_header=1;
else
    warning('function ft_read_header() not found: download FieldTrip and/or adjust path');
end
if 2==exist('sopen');
    has_sopen=1;
else
    warning('function sopen() not found: download BIOSIG and/or adjust path');
end
    
    
if has_ft_read_header
    isftReadable=0;
    try
        % readable by FieldTrip?
        h=ft_read_header(fname);
        isftReadable=1;
    catch
        ; % can't read
    end
end
if ~isftReadable & has_sopen
    isBiosigReadable=0;
    try
        % readable by biosig?
        h=sopen(fname);
        isBiosigReadable=1;
        sclose(h);
    catch
        ; % can't read
    end
end
    
if isftReadable
    if VERBOSE; disp('read with FieldTrip'); end
    h=ft_read_header(fname);    
    p.header=h;
    p.read_with='FieldTrip';
    p.sr=h.Fs;
    if nargout>1;
        data=ft_read_data(fname)';
    end
elseif isBiosigReadable
    if VERBOSE; disp('read with Biosig'); end
    h=sopen(fname);
    p.header=h;
    p.read_with='BIOSIG';
    p.sr=h.SampleRate;
    if nargout>1;
        data=sread(h)';
    end
    sclose(h);
else
    ismatfile=0;
    try
        % .mat file?
        S=whos('-file',data);
        if numel(S)>1
            if nargout==2
                [p,data]=nt_read_data([fname,'.mat']);
            else
                [p,data]=nt_read_data([fname,'.mat']);
            end
        end
    catch
        disp(['File >',fname,'< is not a matlab file, and FieldTrip and BIOSIG can''t read it']);
    end
end
    

