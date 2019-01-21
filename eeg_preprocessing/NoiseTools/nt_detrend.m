function [y,w,r,regressline]=nt_detrend(x,order,w,basis,thresh,niter)
%[y,w,r]=nt_detrend(x,order,w,basis,thresh,niter) - robustly remove trend
% 
%  y: detrended data
%  w: updated weights
%  r: basis matrix used
%
%  x: raw data
%  order: order of polynomial or number of sin/cosine pairs
%  w: weights
%  basis: 'polynomials' [default] or 'sinusoids', or user-provided matrix
%  thresh: threshold for outliers [default: 3 sd]
%  niter: number of iterations [default: 5]
%
% Noise tools
% See nt_regw().
%
% The data are fit to the basis using weighted least squares. The weight is
% updated by setting samples for which the residual is greater than 'thresh' 
% times its std to zero, and the fit is repeated at most 'niter'-1 times.
%
% The choice of order (and basis) determines what complexity of the trend
% that can be removed.  It may be useful to first detrend with a low order
% to avoid fitting outliers, and then increase the order.
%
% Examples:
% Fit linear trend, ignoring samples > 3*sd from it, and remove:
%   y=nt_detrend(x,1); 
% Fit/remove polynomial order=5 with initial weighting w, threshold = 4*sd:
%   y=nt_detrend(x,5,w,[],4);
% Fit/remove linear then 3rd order polynomial:
%   [y,w]=nt_detrend(x,1);
%   [yy,ww]=nt_detrend(y,3);
%
nt_greetings;

%% arguments
if nargin<2; error('!'); end
if nargin<3; w=[]; end
if nargin<4||isempty(basis); basis='polynomials'; end
if nargin<5||isempty(thresh); thresh=3; end
if nargin<6||isempty(niter); niter=4; end

if thresh==0; error('thresh=0 is not what you want...'); end % common mistake

dims=size(x);
x=x(:,:); % concatenates dims >= 2
w=w(:,:);

if size(w,2)==1; w=repmat(w,1,size(x,2)); end

%% regressors
if isnumeric(basis)
    r=basis;
else
    switch basis
        case 'polynomials'
            r=zeros(size(x,1),numel(order));
            lin=linspace(-1,1,size(x,1));
            for k=1:order
                r(:,k)=lin.^k;
            end
        case 'sinusoids'
            r=zeros(size(x,1),numel(order)*2);
            lin=linspace(-1,1,size(x,1));
            for k=1:order
                r(:,2*k-1)=sin(2*pi*k*lin/2);
                r(:,2*k)=cos(2*pi*k*lin/2);
            end
        otherwise
            error('!');
    end
end
%r=nt_normcol(nt_demean(r));

%% remove trends

% The tricky bit is to ensure that weighted means are removed before
% calculating the regression (see nt_regw).

for iIter=1:niter
    
    %disp(iIter); 
    %nt_whoss;
    
    % weighted regression on basis
    [~,y]=nt_regw(x,r,w);
    
    % find outliers
    d=x-y; 
    if ~isempty(w); d=bsxfun(@times,d,w); end
    ww=ones(size(x));
    ww(find(abs(d)>thresh*repmat(std(d),size(x,1),1))) = 0;
    clear d
    
    % update weights
    if isempty(w); 
        w=ww;
    else
        w=min(w,ww);
    end
    clear ww;
    
end

regressline = reshape(y,dims);

%y=x-y;
y=reshape(y,dims);
w=reshape(w,dims);

%new
x=reshape(x,dims);
y=x-y;

if ~nargout
    % don't return, just plot
    figure(1); clf;
    subplot 411; plot(x); title('raw');
    subplot 412; plot(y); title('detrended');
    subplot 413; plot(x-y); title('trend');
    subplot 414; nt_imagescc(w'); title('weight');
    clear y w r
end





%% test code
if 0
    % basic
    x=(1:100)'; x=x+ randn(size(x));
    y=nt_detrend(x,1);
    figure(1); clf; plot([x,y]);
end
if 0
    % detrend biased random walk
    x=cumsum(randn(1000,1)+0.1);
    y=nt_detrend(x,3,[]);
    figure(1); clf; plot([x,y]); legend('before', 'after');
end
if 0
    % weights
    x=linspace(0,100,1000)';
    x=x+3*randn(size(x));
    x(1:100,:)=100;
    w=ones(size(x)); w(1:100,:)=0;
    y=nt_detrend(x,3,[],[],100);
    yy=nt_detrend(x,3,w);
    figure(1); clf; plot([x,y,yy]); legend('before', 'unweighted','weighted');
end
if 0
    [p,x]=nt_read_data('/data/meg/theoldmanandthesea/eeg/mc/MC_aespa_speech_45.mat'); x=x'; x=x(:,1:128); %x=x(1:10000,:);
    %[p,x]=nt_read_data('/data/meg/arzounian/ADC_DA_140521_p20/ADC_DA_140521_p20_01_calib'); x=x'; x=x(1:10000,:);
    
    x=nt_demean(x);
    nt_detrend(x,10);   
end

