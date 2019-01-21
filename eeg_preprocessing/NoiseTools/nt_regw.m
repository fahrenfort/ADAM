function [b,z]=nt_regw(y,x,w)
%[b,z]=nt_regw(y,x,w) - weighted regression
%
%  b: regression matrix (apply to x to approximate y)
%  z: regression (x*b)
%
%  y: data
%  x: regressor
%  w: weight to apply to y
%
%  w is either a matrix of same size as y, or a column vector to be applied
%  to each column of y
%
% NoiseTools

PCA_THRESH=0.0000001; % discard dimensions of x with eigenvalue lower than this

if nargin<3; w=[]; end
if nargin<2; error('!'); end

%% check/fix sizes
m=size(y,1);
x=nt_unfold(x);
y=nt_unfold(y);
if size(x,1)~=size(y,1); 
    disp(size(x)); disp(size(y)); error('!'); 
end

%% save weighted mean
if nargout>1
    mn=y-nt_demean(y,w);
end

%%
if isempty(w) 
    %% simple regression
    xx=nt_demean(x);
    yy=nt_demean(y);
    [V,D]=eig(xx'*xx); V=real(V); D=real(D);
    topcs=V(:,find(D/max(D) > PCA_THRESH)); % discard weak dims
    xxx=xx*topcs;
    b=(yy'*xxx) / (xxx'*xxx); b=b';
    if nargout>1; z=nt_demean(x,w)*topcs*b; z=z+mn; end
else
    %% weighted regression
    if size(w,1)~=size(x,1); error('!'); end
    if size(w,2)==1; 
        %% same weight for all channels
        if sum(w(:))==0; 
            %warning('weights all zero');
            b=0;
        else
            yy=nt_demean(y,w).*repmat(w,1,size(y,2)); 
            xx=nt_demean(x,w).*repmat(w,1,size(x,2));  
            [V,D]=eig(xx'*xx); V=real(V); D=real(D); D=diag(D);
            topcs=V(:,find(D/max(D) > PCA_THRESH)); % discard weak dims
            xxx=xx*topcs;
            b=(yy'*xxx) / (xxx'*xxx); b=b';
        end
        if nargout>1; z=nt_demean(x,w)*topcs*b; z=z+mn; end
    else
        %% each channel has own weight
        if size(w,2) ~= size(y,2); error('!'); end 
        if nargout; z=zeros(size(y)); end
        for iChan=1:size(y,2)
            if sum(w(:,iChan))==0; %disp(iChan); 
                %warning('weights all zero'); 
                b=zeros(size(y,2));
            else
                yy=nt_demean(y(:,iChan),w(:,iChan)) .* w(:,iChan); 
                x=nt_demean(x,w(:,iChan)); % remove channel-specific-weighted mean from regressor
                xx=x.*repmat(w(:,iChan),1,size(x,2)); 
                [V,D]=eig(xx'*xx); V=real(V); D=real(diag(D));
                topcs=V(:,find(D/max(D) > PCA_THRESH)); % discard weak dims
                xxx=xx*topcs;
                b(iChan,1:size(topcs,2))=(yy'*xxx) / (xxx'*xxx); 
            end
            if nargout>1; z(:,iChan)=x*(topcs*b(iChan,1:size(topcs,2))') + mn(:,iChan); end
        end
    end             
end

%%
if nargout>1;
    z=nt_fold(z,m);
end

%% test code
if 0
    % basic
    x=randn(100,10); y=randn(100,10); 
    b1=nt_regw(x,y); b2=nt_regw(x,x); b3=nt_regw(x,y,ones(size(x))); 
    figure(1); subplot 131; nt_imagescc(b1); subplot 132; nt_imagescc(b2); subplot 133; nt_imagescc(b3);
end
if 0
    % fit random walk
    y=cumsum(randn(1000,1)); x=(1:1000)'; x=[x,x.^2,x.^3];
    [b,z]=nt_regw(y,x); 
    figure(1); clf; plot([y,z]);
end
if 0
    % weights, random
    y=cumsum(randn(1000,1)); x=(1:1000)'; x=[x,x.^2,x.^3];
    w=rand(size(y));
    [b,z]=nt_regw(y,x,w); 
    figure(1); clf; plot([y,z]);
end
if 0
    % weights, 1st vs 2nd half
    y=cumsum(randn(1000,1))+1000; x=(1:1000)'; x=[x,x.^2,x.^3];
    w=ones(size(y)); w(1:500,:)=0;
    [b,z]=nt_regw(y,x,w); 
    figure(1); clf; plot([y,z]);
end
if 0
    % multichannel
    y=cumsum(randn(1000,2)); x=(1:1000)'; x=[x,x.^2,x.^3];
    w=ones(size(y)); 
    [b,z]=nt_regw(y,x,w); 
    figure(1); clf; plot([y,z]);
end

