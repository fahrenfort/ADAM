function y=nt_video_sns(x,nneighbors)
%y=nt_video_sns(x,nneighbors) - apply SNS locally
%
%  y: processed data
%
%  x: data to process (time * nrows * ncols)
%  nneighbors: number of neighbors to incluce for each pixel
%
% Each channel is projected on its nneighbors closest neighbors, 
% and the channel is replaced by its projection.

if nargin<1; error('!'); end

if nargin<2||isempty(nneighbors);
    disp('default nneighbors = 10');
    nneighbors=10;
end

[nframes,nrows,ncols]=size(x);
x=x(:,:); 
y=nan(size(x));

closest=nt_proximity([nrows,ncols],nneighbors);
for iPixel=1:size(closest,1)
    %disp(iPixel)
    xx=x(:,closest(iPixel,:));
    %[~,a]=nt_regw(x(:,iPixel),xx);
    %y(:,iPixel)=a;
    [V,D]=eig(xx'*xx); V=real(V); D=real(D);
    PCA_THRESH=10^-8;
    topcs=V(:,find(D/max(D) > PCA_THRESH)); % discard weak dims
    xxx=xx*topcs;
    b=( x(:,iPixel)'*xxx ) / (xxx'*xxx);    
    y(:,iPixel)=xxx*b';
end
y=reshape(y,[nframes,nrows,ncols]);


