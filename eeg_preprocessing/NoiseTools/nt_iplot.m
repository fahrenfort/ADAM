function [hh,ii]=nt_iplot(name)
%[hh,ii]=nt_iplot(fname) - plot data file based on index
%
%  hh, ii: header info and index
%
%  name: file name
%
%  nt_iplot(name): just plot the data, mean removed
%
% See nt_index.
%
% NoiseTools
nt_greetings;

if nargin<1 || isempty(name);
    error('!');
end

% check 'name'
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
    error('...is neither file nor directory');
end

if 2~=exist(name)  &&  ~strcmp('.ds',EXT);
    disp('can only handle files');
    disp(['  >',name,'<']);
end

% index directory 
idxDir=[PATHSTR,filesep,'nt_idx'];
if 7 ~= exist(idxDir); 
    nt_index(name);
end

% index file
idxName=[idxDir,filesep,NAME,EXT,'.idx'];
hhh.idxName=idxName;
if ~2==exist(idxName); 
    nt_index(name); 
end

load ('-mat',idxName);  % loads hh, ii

%figure(1); clf
t=0:size(ii.min,1)-1;
if ~isempty(hh.sr) && isfield (ii,'scale'); 
    t=t*ii.scale/hh.sr;
    xlab=['seconds (',num2str(ii.nsamples),' samples)'];
else
    xlab='samples';
end
plot(t,nt_demean([ii.min,ii.max]));
xlabel(xlab);
xlim([t(1) t(end)]);
title(NAME,'interpreter','none');

if ~nargout
    hh=[];ii=[];
end

