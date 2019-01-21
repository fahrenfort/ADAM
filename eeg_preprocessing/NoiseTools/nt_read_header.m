function [h,readwith]=nt_read_header(fname)
%[h,readwith]=nt_read_header(fname,flag) - read data from file
%
%  h: header
%  readwith: fieldtrip or biosemi
% 
%  fname: file to read
%  
% 
if nargin < 1 ; error('!'); end
if ~isa(fname, 'char'); 
    error('filename is not char string');
end

[PATHSTR,NAME,EXT] = fileparts(fname);
if isempty(PATHSTR); 
    fname=[pwd,filesep,fname]; % safe to use exist
end

if 7==exist(fname);
    disp(fname);
    error('...is directory!')
end

if exist(fname,'file')~=2; 
    disp(fname)
    error('...not found');
end

if numel(fname)>4 & fname(end-3:end)=='.mat' % intercept matlab files
    disp(fname)
    error('...is mat file');
end   
         
% select file reader among those available
persistent nt_read_header_readwith
if isempty(nt_read_header_readwith)
    if 2==exist('ft_read_header');
        nt_read_header_readwith='fieldtrip';
    else
        warning('function ft_read_header() not found: download FieldTrip and/or adjust path');
    end
end
if isempty(nt_read_header_readwith)
    if 2==exist('sopen');
        nt_read_header_readwith='biosig';
    else
        warning('function sopen() not found: download BIOSIG and/or adjust path');
    end
end
if isempty(nt_read_header_readwith); error('no reading functions available'); end

% read the header
if strcmp(nt_read_header_readwith,'fieldtrip');
    try
        h=ft_read_header(fname);
    catch
        h=[];
    end
else
    try
        h=sopen(fname);
        sclose(h);
    catch
        h=[];
    end
end

readwith=nt_read_header_readwith;

