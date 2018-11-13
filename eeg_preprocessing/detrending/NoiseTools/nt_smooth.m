function x=nt_smooth(x,T,nIterations,nodelayflag)
%y=nt_smooth(x,T,nIterations,nodelayflag) - smooth by convolution with square window
%
%  y: smoothed data
% 
%  x: data to smooth
%  T: samples, size of window (can be fractionary)
%  nIterations: number of iterations of smoothing operation (large --> gaussian kernel)
%  nodelayflag: if true, compensate for delay [default:false]
%
nt_greetings;

if nargin<4||isempty(nodelayflag); nodelayflag=0; end
if nargin<3||isempty(nIterations); nIterations=1; end
if nargin<2; help nt_smooth ; error; end

if ndims(x)>4; error('!'); end

integ=floor(T);
frac=T-integ;

if integ>=size(x,1);
    x=repmat(mean(x),[size(x,1),1,1,1]);
    return;
end

% remove onset step
mn=mean(x(1:(integ+1),:,:),1);
x=bsxfun(@minus,x,mn);

if nIterations==1 && frac==0;
    % faster
    x=cumsum(x);
    x(T+1:end,:)=x(T+1:end,:)-x(1:end-T,:);
    x=x/T;
else
    % filter kernel
    B=[ones(integ,1);frac]/T;
    for k=1:nIterations-1
        B=conv(B,[ones(integ,1);frac]/T);
    end
    x=filter(B,1,x);    
end

if nodelayflag
    shift=round(T/2*nIterations); %[shift n*T]
    x=[x(shift+1:end,:,:,:); zeros(shift,size(x,2),size(x,3),size(x,4))];
end


% restore DC
x=bsxfun(@plus,x,mn);