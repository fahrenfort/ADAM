 function x=nt_dss_repeat_cluster(x,nkeep,cluster_size)
% y=nt_dss_repeat_cluster(x,nkeep,cluster_size) - dss within clusters
% to emphasize repeatability
%
%   y: denoised matrix
%
%   x: matrix  to denoise (time * channels * trials)
%   nkeep: number of components to keep for each cluster
%   cluster_size: target cluster size [default: size(x,1)/2]
%
% NoiseTools
nt_greetings;


if nargin<3; cluster_size=[]; end
if nargin<2; nkeep=ceil(cluster_size/2); end

if ndims(x)~=3; error('!'); end

x0=x;

if isempty(cluster_size); cluster_size=round(size(x,1)/2); end

if size(x,2)<=cluster_size
    todss=nt_dss1(x);
    fromdss=pinv(todss);
    nkeep=min(nkeep,size(todss,2));
    x=nt_mmat(x,(todss(:,1:nkeep)*fromdss(1:nkeep,:)));
else
    NCLUSTERS=2;
    [C,A]=vl_kmeans(nt_unfold(x),NCLUSTERS,'algorithm', 'elkan','initialization','plusplus','numrepetitions', 100);
    if numel(find(A==1)) && numel(find(A==2))
        xA=nt_dss_repeat_cluster(x(:,find(A==1),:),nkeep,cluster_size);
        xB=nt_dss_repeat_cluster(x(:,find(A==2),:),nkeep,cluster_size);
        x(:,find(A==1),:)=xA;
        x(:,find(A==2),:)=xB;
    end % else no split, return
end
x1=x;


verbose=0;
if nargout==0 || verbose;
    disp(['cluster size: ', num2str(size(x,2)), ',  power ratio: ',num2str(nt_wpwr(x1)/nt_wpwr(x0))]);
end


