 function x=nt_sns_cluster(x,nneighbors,cluster_size)
% y=nt_sns_cluster(x,nneigbors,cluster_size) - sensor noise suppression within clusters
%
%   y: denoised matrix
%
%   x: matrix  to denoise (time * channels * repeats * ...)
%   nneighbors: number of channels to use in projection
%   cluster_size: target cluster size [default: size(x,1)/2]
%
% NoiseTools
nt_greetings;


if nargin<3; cluster_size=[]; end
if nargin<2; error('!'); end

sz=size(x);
x=x(:,:,:); % merge higher dimensions into repeats
x=nt_unfold(x);
x0=x;

if isempty(cluster_size); cluster_size=round(size(x,1)/2); end

if size(x,2)<=cluster_size
    x=nt_sns(x,nneighbors);
else
    NCLUSTERS=2;
    [C,A]=vl_kmeans(x,NCLUSTERS,'algorithm', 'elkan','initialization','plusplus','numrepetitions', 100);
    if numel(find(A==1)) && numel(find(A==2))
        xA=nt_sns_cluster(x(:,find(A==1)),nneighbors,cluster_size);
        xB=nt_sns_cluster(x(:,find(A==2)),nneighbors,cluster_size);
        x(:,find(A==1))=xA;
        x(:,find(A==2))=xB;
    end % else no split, return
end
x1=x;

x=nt_fold(x,sz(1));
x=reshape(x,sz);

verbose=0;
if nargout==0 || verbose;
    disp(['cluster size: ', num2str(sz(2)), ',  power ratio: ',num2str(nt_wpwr(x1)/nt_wpwr(x0))]);
end


