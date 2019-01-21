function [IDX,TODSS,SCORE,COVS]=nt_cluster_jd(x,dsr,smooth,flags,init,verbose, depth,N)
%[IDX,todss,SCORE,COVS]=nt_cluster_jd(x,dsr,smooth,flags,init,verbose) - cluster with joint diagonalization
%
%  IDX: cluster ownership (IDX{1}: low amp, IDX{2{: high amp)
%  TODSS: DSS matrix (1st column --> discriminating component)
%  SCORE: score (smaller means better contrast)
%  COVS: covariance for each cluster
%
%  x: data (time*channel% s)
%  dsr: downsample ratio for cross product series
%  smooth: further smoothing of cross-product series
%  flags: see below
%  init: provide initial clustering
%  verbose: display & plot (default=no)
%  depth: cluster recursively into 2^depth clusters
%  N: target number of clusters [default: depth^2]
%
% Flags:
%  'norm', 'norm2': give each slice the same weight
%  'amp', 'pwr': cluster amplitude or power instead of log (default)
% See nt_bias_cluster, nt_cluster1D


if nargin<2; error('!'); end
if nargin<3 ||isempty(smooth); smooth=1; end
if nargin<4 ||isempty(flags); flags=[]; end
if nargin<5; init=[]; end
if nargin<6||isempty(verbose); verbose=0; end
if nargin<7||isempty(depth); depth=1; end
if nargin<8||isempty(N); N=2^depth; end

if ndims(x)>2 || size(x,2) ==1;
    error('should be 2D matrix');
end

if depth>1;
    % split into 2 clusters
    I=nt_cluster_jd(x,dsr,smooth,flags,init,verbose,1);
    
    % recurse on first
    if numel(I{1})>2*dsr; 
        I1=nt_cluster_jd(x(I{1},:),dsr,smooth,flags,init,verbose,depth-1); % recurse
    else
        I1={(1:numel(I{1}))}; % too small
    end
    
    % recurse on second (if exists)
    if numel(I)>1; 
        if numel(I{2})>2*dsr; I2=nt_cluster_jd(x(I{2},:),dsr,smooth,flags,init,verbose,depth-1); else I2={(1:numel(I{2}))}; end
    end
    
    % resolve the cluster indices
    IDX1={};
    for k=1:numel(I1)
        IDX1=[IDX1, I{1}(I1{k})];
    end
    IDX2={};
    if numel(I)>1
        for k=1:numel(I2)
            IDX2=[IDX2, I{2}(I2{k})];
        end
    end
    IDX=[IDX1 IDX2];
    checkindex(IDX,size(x,1))
    
    while numel(IDX)>N;
        % merge clusters
        COVS=[];
        for k=1:numel(IDX)
            COVS{k}=nt_cov(x(IDX{k},:))/size(x(IDX{k},:),1);
        end
        B=covdists(COVS);
        [a,idx]=min(B(:));
        [k1,k2]=ind2sub([size(B,1), size(B,1)],idx);
        
        %figure(1); clf; nt_imagescc(B);title (num2str([max(B(:)), a, k1, k2, size(B,1)])); pause
        IDX{k1}=[IDX{k1};IDX{k2}]; 
        IDX(k2)=[];
        checkindex(IDX,size(x,1))
    end
    if nargout>1;
        c0=nt_cov(x)/size(x,1);
        for k=1:numel(IDX)
            c1=nt_cov(x(IDX{k},:))/size(x(IDX{k},:),1);
            [TODSS{k},pwr0,pwr1]=nt_dss0(c0,c1);
            SCORE(k,1:numel(pwr1))=pwr1./pwr0;
            COVS{k}=c1;
        end
    end
    return
end    

%{
 Calculate the time series of cross products (terms of the covariance matrix).
 This time series has coarser temporal resolution than x by a factor dsr.
%}
[xx,ind]=nt_xprod(x,'lower',dsr);
if 0
    disp([num2str(size(xx,2)), ' crossproducts']);
    nt_whoss;
end

% figure(2); clf;
% subplot 211;
% plot(xx)

% option: give each slice the same weight (counters amplitude variations)
if find(strcmp(flags,'norm'))
    xx=nt_normrow(xx);
end
if find(strcmp(flags,'norm2'))
    xx=norm2(xx,size(x,2),ind);
end

% subplot 212; 
% plot(xx); 
% pause;

xx=nt_smooth(xx,smooth,[],1);

%{
Cluster each column the time series of cross products, 
choose the column with best score (reduction in energy), 
and use it's cluster index to initialize the first JD analysis.
%}

% initial clustering, DSS
if isempty(init)
    [C,A,score]=nt_cluster1D(xx); % cluster all columns of cross products
    [~,idx]=min(score); % select column with best score (tightest clusters)
    A=A(:,idx); 
        
    % upsample the cluster ownership index so we can apply it to x
    A=repmat(A',[dsr,1]);
    A=A(:);
    A(end:size(x,1))=A(end);
    IDX{1}=find(A==0);
else
    IDX{1}=init;
end

if isempty(IDX{1}) % clustering failed, return just one cluster
    IDX{1}=1:size(x,1);
    TODSS{1}=nan;
    SCORE{1}=nan;
    COVS{1}=nt_cov(x);
    return
end
   
c0=nt_cov(x);
c1=nt_cov(x(IDX{1},:));
[todss,pwr0,pwr1]=nt_dss0(c0,c1);
z=nt_mmat(x,todss(:,[1 end])); % keep only first and last components 

PLOT_FIG2=0;
if PLOT_FIG2
    figure(2);  clf; set(gcf, 'name','in nt_cluster_jd');
    A=zeros(size(x,1),1); A(IDX{1})=1;
    subplot 511; plot(x); title('data');
    subplot 512; plot(A,'.-'); title('initial cluster map');
    subplot 513; plot(z(:,1)); title('initial DSS1');
    subplot 514; plot(z(:,2)); title('initial DSS2');
    drawnow; pause;
end

% iterate until stable
old_IDX=IDX{1};
for k=1:10

    [zz,ind]=nt_xprod(z,[],dsr);
    zz=zz(:,1:2);       % keep only the squares
    
    if find(strcmp(flags,'pwr')); % cluster in power 
        [C,A]=nt_cluster1D(zz);
        [~,idx]= max(abs(diff(log2(C+eps)))); % choose first or last
    elseif find(strcmp(flags,'amp')); % cluster in amplitude 
        [C,A]=nt_cluster1D(sqrt(zz));
        [~,idx]= max(abs(diff(log2(C+eps)))); % choose first or last
    else  % cluster in log domain
        [C,A]=nt_cluster1D(log2(zz+eps));
        [~,idx]= max(abs(diff(C))); % choose first or last
    end
    A=A(:,idx);
    %disp(C);
    C=C(:,idx);
    %disp(C); pause
    if C(1)<C(2); A=1-A; end % ensure that first cluster has low amplitude
    
    A=double(nt_smooth(A,smooth, [],1)>=1/smooth); % extend ownership to include effect of smoothing

    % upsample the cluster ownership index so we can apply it to x
    A=repmat(A',[dsr,1]); % upsample 
    A=A(:); 
    A(end:size(x,1))=A(end);
    IDX{1}=find(A==0); % 0: low values, 1: high values
    
    if isempty(IDX{1}) % clustering failed, return just one cluster
        IDX{1}=1:size(x,1);
        TODSS{1}=nan;
        SCORE{1}=nan;
        COVS{1}=nt_cov(x);
        return
    end

    % DSS to contrast clusters
    c0=nt_cov(x)/size(x,1);
    c1=nt_cov(x(IDX{1},:))/size(x(IDX{1},:),1);
    [todss,pwr0,pwr1]=nt_dss0(c0,c1);
    z=nt_mmat(x,todss(:,[1 end])); % keep first and last

    if ~nargout||verbose; 
        disp(['low amp cluster: ', num2str((100*numel(IDX{1})/size(x,1)), 2), ' % of samples, power ratio: ' num2str(pwr1(end)/pwr0(end), 3)]); 
        disp(['hi amp cluster: ', num2str((100-100*numel(IDX{1})/size(x,1)), 2), ' % of samples, power ratio: ' num2str(pwr1(1)/pwr0(1), 3)]); 
    end

    if PLOT_FIG2
        figure(2);  
        subplot 515; plot(A,'.-'); title('final cluster map'); pause
    end
    if all(size(old_IDX)==size(IDX{1})) && all(old_IDX==IDX{1}); break; end
    old_IDX=IDX{1};
end 
IDX{2}=setdiff((1:size(x,1))', IDX{1});


% final DSS
c0=nt_cov(x)/size(x,1);
c1=nt_cov(x(IDX{1},:))/size(x(IDX{1},:),1);
COVS{1}=c1;
[TODSS{1},pwr0,pwr1]=nt_dss0(c0,c1);
SCORE(1,1:numel(pwr1))=pwr1./pwr0;
c1=nt_cov(x(IDX{2},:))/size(x(IDX{2},:),1);
COVS{2}=c1;
[TODSS{2},pwr0,pwr1]=nt_dss0(c0,c1);
SCORE(2,1:numel(pwr1))=pwr1./pwr0;

if nargout==0||verbose;
    
    % no output, just plot

    z1=nt_mmat(x,TODSS{1}(:,1));

    figure(101); clf ;
    subplot 221;
    plot(pwr1./pwr0,'.-'); xlabel('component'); ylabel('score'); title('DSS cluster vs all');
    subplot 222;
    wsize=min(1024,size(z1,1));
    hold on
    nt_spect_plot(z1/sqrt(mean(z1(:).^2)),wsize,[],[],1);
    nt_spect_plot(x/sqrt(mean(x(:).^2)),wsize,[],[],1);
    xlim([0 .5]);
    nt_linecolors([],[1 3 2]);
    legend('cluster','all'); legend boxoff
    hold off

    z=nt_mmat(x,todss); 
    z=nt_normcol(z);
    subplot 223; imagescc(nt_cov(z(IDX{1},:))); title('cluster 1'); 
    subplot 224; imagescc(nt_cov(z)-nt_cov(z(IDX{1},:))); title('cluster 2');

    
    figure(102); clf
    if 0
        subplot 211;
        plot(x); hold on
        xx=x; xx(IDX{1},:)=nan;
        plot(xx,'k');
        axis tight
        title('black: cluster [high amp]');
        subplot 212;
        plot(z1); axis tight
        title('first DSS component');
    else
        subplot 311;
        plot(x); hold on
        xx=x; xx(IDX{1},:)=nan;
        plot(xx,'k');
        axis tight
        title('black: cluster [high amp]');
        subplot 312;
        plot(z1); axis tight
        title('DSS 1');
        subplot 313;
        nt_sgram(z1,128,1); axis tight
        title('DSS 1');
    end
    
    if 0 
        figure(105); clf
        nt_sgram(z1,1024,32,[],1);
        title('DSS1');
    end
    if nargout==0; clear IDX SCORE TODSS; end
    
end

% can't rememember what this is supposed to do...
function y=norm2(x,nchans,ind)
[I,J]=ind2sub([nchans,nchans],ind); % linear --> matrix indices
for k=1:size(x,1)
    a=x(k,1:nchans);
    b=sqrt(a(I).*a(J));
    y(k,:)=x(k,:)./b;
end

% matrix of covariance distances
function B=covdists(C) % B: matrix of distances, C: array of covariance matrices
B=nan(numel(C));
CC=zeros(size(C{1}));
for k=1:numel(C); CC=CC+C{k}; end
for k=1:numel(C)
    for j=1:k-1
        [E]=eig(abs(C{j}-C{k}),CC);
        B(j,k)=max(abs(log2(E)));
        B(k,j)=B(j,k);
    end
end
    
function checkindex(IDX,n)
a=zeros(n,1);
for k=1:numel(IDX)
    if any(a(IDX{k})); 
        error('!');
    end
    a(IDX{k})=1;
end
    
    