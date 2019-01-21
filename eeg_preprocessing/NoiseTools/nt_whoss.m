function varargout=nt_whoss
%size=nt_whoss - total Gbytes used by variables
%
%  size: number of Gbytes
%
% If nargout==0, display Gbytes used

s=evalin('caller', 'whos'); % lists variables in caller's workspace

x=cell(1,numel(s));         
[x{:}]=deal(s.bytes);       % transfer to cell array
size=sum(cat(1,x{:}));      % transfer to array and sum

size=size/(2^30);           % bytes --> Gbytes
nfiles=numel(fopen('all')); % number of open files

if nargout == 0
    ST=dbstack;
    if numel(ST)>1;
        disp(['Gbytes used: ', num2str(size), ', open files: ', num2str(nfiles), ...
            ', ', ST(2).file,' line ', num2str(ST(2).line)]);
    else
        disp(['Gbytes used: ', num2str(size), ', open files: ', num2str(nfiles)]);
    end
else
    varargout{1}=size;
    varargout{2}=nfiles;
end

