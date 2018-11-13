function y=nt_sparse_filter(x,T,A)
%y=nt_sparse_filter(x,T,A) - convolve multichannel data with sparse impulse response
%
%  y: result
%
%  x: data to convolve (columnwise)
%  T: times of non-zero samples of IR (can be fractionnary)
%  A: amplitude of non-zero samples of IR (default: all 1)
%
% T and A together describe the impulse response.  A must have same size
% as T.
%
% If T is a column vector it is applied to all columns of x.
%
% If T is a 2D matrix or 1D cell array each column or cell is applied to a
% column of x (dimensions must fit).
%
% If T is a 2D cell array it is interpreted as defining a multichannel 
% impulse response. The cell on row I and column J is applied to column I
% of x and the result is added to column J of y.
%
% Fractionary lags are implemented by linear interpolation.

if nargin<3; A=[]; end
if nargin<2; error('!'); end

% check A or create
%if isvector(T); T=T(:); A=A(:); end
if isnumeric(T)
    if isempty(A); A=ones(size(T)); end
    if ~isnumeric(A); error('!'); end
    if size(A) ~= size(T); error('!'); end
    if size(A,2) ~= 1 || (size(A,2) ~= size(x,2) && size(A,2)~=1); error('!'); end
elseif iscell(T)
    if isempty(A); 
        A=cell(size(T));
        for iCell=1:numel(T); A{iCell}=ones(size(T{iCell})); end
    end
    if ~iscell(A); error('!'); end
    if size(A) ~= size(T); error('!'); end
    if size(A,1)~=size(x,2); error('number of rows of cell array of IRs should equal number of columns of x'); end
end


% handle negative lags
if isnumeric(T)
    a=min(T(:));
    if a<0
        x=[zeros(ceil(-a)),size(x,2)];
        T=T+ceil(-a);
    end
elseif iscell(T)
    a=min(T{1,1}); 
    for iCell=1:numel(T)
        a=min(a,min(T{iCell}));
    end
    if a<0
        x=[zeros(ceil(-a)),size(x,2)];
        for iCell=1:numel(T)
            T{iCell}=T{iCell}+ceil(-a);
        end
    end
else; error('!'); end % neither matrix nor cell

% estimate nrows of result
if isnumeric(T)
    b=max(T(:));
else
    b=max(T{1,1});
    for iCell=1:numel(T)
        b=max(b,max(T{iCell}));
    end
end
b=ceil(b);

[nsamples,ncolsx]=size(x);

if isnumeric(T) && size(T,2)==1 
    % apply same IR to all columns of x
    y=zeros(nsamples+b, ncolsx);
    for iPulse=1:numel(T);
        t=T(iPulse);
        integT=floor(t); 
        fracT=t-integT;
        y(integT+(1:nsamples),:) = y(integT+(1:nsamples),:) + x*(1-fracT)*A(iPulse);
        if fracT; y(1+integT+(1:nsamples),:) = y(1+integT+(1:nsamples),:) + x*fracT*A(iPulse); end
    end
end

if isnumeric(T) && size(T,2)>1
    % apply different IR to each column of x
    y=zeros(nsamples+b, ncolsx);
    for iCol=1:ncolsx
        for iPulse=1:numel(T(:,iCol));
            t=T(iPulse,iCol);
            integT=floor(t); 
            fracT=t-integT;
            1-fracT
            y(integT+(1:nsamples),iCol) = y(integT+(1:nsamples),iCol) + x(:,iCol)*(1-fracT)*A(iPulse,iCol);
            if fracT; y(1+integT+(1:nsamples),iCol) = y(1+integT+(1:nsamples),iCol) + x(:,iCol)*fracT*A(iPulse,iCol); end
        end
    end
end

if iscell(T) && numel(T)==ncolsx
    % apply one IR to each column, no cross-terms
    y=zeros(nsamples+b, ncolsx);
    for iCol=1:ncolsx
        for iPulse=1:numel(T{iCol});
            t=T{iCol}(iPulse);
            integT=floor(t); 
            fracT=t-integT;
            y(integT+(1:nsamples),iCol) = y(integT+(1:nsamples),iCol) + x(:,iCol)*(1-fracT)*A{iCol}(iPulse);
            if fracT; y(1+integT+(1:nsamples),iCol) = y(1+integT+(1:nsamples),iCol) + x(:,iCol)*fracT*A{iCol}(iPulse); end
        end
    end
end

if iscell(T) && size(T,1)==ncolsx 
    ncolsy=size(T,2);
    % full multichannel IR with cross-terms
    y=zeros(nsamples+b, ncolsy);
    for iRow=1:ncolsx
        for iCol=1:ncolsy
            for iPulse=1:numel(T{iRow,iCol})
                t=T{iRow,iCol}(iPulse);
                integT=floor(t); 
                fracT=t-integT;
                y(integT+(1:nsamples),iCol) = y(integT+(1:nsamples),iCol) + x(:,iRow)*(1-fracT)*A{iRow,iCol}(iPulse);
                if fracT; y(1+integT+(1:nsamples),iCol) = y(1+integT+(1:nsamples),iCol) + x(iRow)*fracT*A{iRow,iCol}(iPulse); end
            end
        end
    end
end


return

% test code
x=randn(1000,1);
y=nt_sparse_filter(x,[0]);
figure(1); clf; subplot 211; plot(x); subplot 212; plot(y);

x=randn(1000,1);
y=nt_sparse_filter(x,[0.5]);
figure(1); clf; subplot 211; plot(x); subplot 212; plot(y);

x=randn(1000,1);
y=nt_sparse_filter(x,[10]);
figure(1); clf; subplot 211; plot(x); subplot 212; plot(y);

x=randn(1000,1);
y=nt_sparse_filter(x,[1:10]);
figure(1); clf; subplot 211; plot(x); subplot 212; plot(y);

x=randn(1000,1);
y=nt_sparse_filter(x,[1:10], ones(1,10));
figure(1); clf; subplot 211; plot(x); subplot 212; plot(y);

x=randn(1000,2);
y=nt_sparse_filter(x,[1:10], ones(1,10));
figure(1); clf; subplot 211; plot(x); subplot 212; plot(y);

x=randn(1000,2);
y=nt_sparse_filter(x,[1:10], ones(2,10));
figure(1); clf; subplot 211; plot(x); subplot 212; plot(y);

x=randn(1000,2);
y=nt_sparse_filter(x,{ones(1,10),ones(1,10)});
figure(1); clf; subplot 211; plot(x); subplot 212; plot(y);

x=randn(1000,1);
y=nt_sparse_filter(x,{ones(1,10),ones(1,10)});
figure(1); clf; subplot 211; plot(x); subplot 212; plot(y);

            
    