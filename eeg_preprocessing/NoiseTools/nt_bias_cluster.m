function [c0,c1,A,todss,pwr0,pwr1]=nt_bias_cluster(x,dsr,flags)
%[c0,c1,A,todss,pwr0,pwr1]=nt_bias_cluster(x,dsr,flags) - cluster covariance
%
%  c0,c1: covariance matrices of clusters
%  A: map of cluster ownership
%  todss,pwr0,pwr1: result of DSS
%
%  x: data (time*channels)
%  dsr: downsample ratio for cross product series
%  flags: 'norm': give each dsr-sized slice the same weight
%
%  See: nt_cluster1D, nt_cluster_jd.

SMOOTH=1;%2; % smooth the time series of cross products

if nargin<3 ||isempty(flags); flags=[]; end
if nargin<2; error('!'); end

if ~exist('vl_kmeans');
    disp('vl_kmeans() not found, download from http://www.vlfeat.org');
end
if ndims(x)>2; 
    error('x should be time*channels');
end

% time series of cross-products
if find(strcmp(flags,'nodiag'))
    [xx,ind]=nt_xprod(x,'nodiag',dsr);
else
    [xx,ind]=nt_xprod(x,'lower',dsr);
end

% smooth
xx=filter(ones(SMOOTH,1),1,xx); 

% give each slice the same weight (counters amplitude variations)
if find(strcmp(flags,'norm'))
    xx=nt_normrow(xx);
end
if find(strcmp(flags,'norm2'))
    xx=norm2(xx,size(x,2),ind);
end

% cluster the time series (2 clusters)
NCLUSTERS=2;
[C,A]=vl_kmeans(xx',NCLUSTERS,'algorithm', 'elkan','initialization','plusplus',...
    'numrepetitions', 100);

% make sure the first cluster is biggest
if numel(find(A==1))<numel(A)/2;
    C=fliplr(C);
    A=3-A;
end

% upsample the cluster ownership index 
A=repmat(A,[dsr,1]);
A=A(:);
A(end:size(x,1))=A(end);

if 1
c0=nt_cov(x(find(A==1),:));
c1=nt_cov(x(find(A==2),:));
else
% full covariance matrices from lower diagonal values
c0=squeeze(nt_lower_to_full(C(:,1)',ind));   
c1=squeeze(nt_lower_to_full(C(:,2)',ind));   
end

% DSS to find components maximally different between clusters
[todss,pwr0,pwr1]=nt_dss0(c0+c1,c1);


if nargout==0;
    % no output, just plot
    disp(['cluster1: ',num2str(100*numel(find(A==1))/numel(A)), '%']);

    z1=nt_mmat(x,todss(:,1));
    z2=nt_mmat(x,todss(:,end));
    z=nt_mmat(x,todss); 
    z=nt_normcol(z);
    e1=mean(z(find(A==1),:).^2);
    e2=mean(z(find(A==2),:).^2);

    figure(100); clf
    plot(x); hold on
    x(find(A==2),:)=nan;
    plot(x,'k');
    axis tight
    title('black: cluster2');
    
    figure(101); clf
    subplot 121;
    plot(pwr1./pwr0,'.-'); xlabel('component'); ylabel('score'); title('DSS cluster A vs B');
    subplot 122;
    nt_spect_plot(z1,1024,[],[],1);
    hold on
    nt_spect_plot(z2,1024,[],[],1);
    xlim([0 .5])
    nt_linecolors([],[1 3]);
    legend('first','last'); legend boxoff
    hold off

    
    figure(102); clf
    subplot 211;
    plot(z1); axis tight
    title('first DSS component')
    subplot 212;
    plot(z2); axis tight
    title('last DSS component');
    
    figure(103); clf
    plot([e1',e2'], '.-'); legend('cluster A', 'cluster B'); title ('power per component');
    
    figure(104);
    subplot 121; imagescc(c0); title('cluster A'); 
    subplot 122; imagescc(c1); title('cluster B');
    
    if 0 
        figure(105); clf
        subplot 211;
        nt_sgram(z1,1024,32,[],1);
        title('first');
        subplot 212;
        nt_sgram(z2,1024,32,[],1);
        title('last');
    end
    clear c0 c1 A todss pwr0 pwr1
    
end

function y=norm2(x,n,ind)
[I,J]=ind2sub([n,n],ind);
for k=1:size(x,1)
    a=x(k,1:n);
    b=sqrt(a(I).*a(J));
    y(k,:)=x(k,:)./b;
end

    
    
    