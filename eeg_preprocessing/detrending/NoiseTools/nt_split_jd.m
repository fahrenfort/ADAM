function [idx,score_vector,todss]=nt_split_jd(x,thresh,depth);
%[idx,score_vector,todss]=nt_split_dss(x,thresh,depth) - segmentation based on joint diagonalization
%
%  idx: index at which to split
%  score_vector: from nt_split
%  todss: DSS matrix
%
%  x: data
%  thresh: truncation threshold for PCA
%  depth: recursion depth

if nargin<3||isempty(depth); depth=1; end
if nargin<2||isempty(thresh); thresh=0; end
if isempty(x); error('!'); end
if ndims(x)>2; error('x should be 2D'); end

[m,n]=size(x);

% initial PCA to remove below-threshold dimensions
topcs=nt_pca0(x);
z=nt_mmat(x,topcs);
keep=find(mean(z.^2)/mean(z(:,1).^2)>thresh);
z=z(:,keep);

if m<=2; warning('m==',num2str(m)); end

idx=ceil(m/2); % initial split into two arbitrary intervals
% iterate until stable
for k=1:10
    c1=nt_cov(z(1:idx,:));
    c0=c1+nt_cov(z(idx+1:end,:));
    todss=nt_dss0(c0,c1);
    zz=nt_mmat(z,todss(:,[1,end]));
    zz=nt_normcol(zz);
    old_idx=idx;
    [idx,score_vector]=nt_split(nt_normcol(nt_demean(zz.^2)));
    %figure(1); clf; subplot 211; plot(zz(:,1));subplot 212; plot(zz(:,2)); idx, pause
    if idx==old_idx; break; end
    disp(num2str([idx,old_idx]));
end
todss=topcs(:,keep)*todss; % to return

if depth>1
    [a]=nt_split_jd(x(1:idx,:),thresh, depth-1);
    [b]=nt_split_jd(x(idx+1:end,:), thresh,depth-1);
    idx=[a,idx,idx+b];
    idx=unique(idx);
end

% disp(['depth, ndims: ' num2str([depth, size(z,2)])])
% disp(['idx: ' num2str([idx])])
% nt_split(nt_normcol(nt_demean(zz.^2)));

disp(['nt_split_jd nargout: ', num2str(nargout)])

if nargout==0;
    disp(['split at ', num2str(idx)]);
    disp(['(%: ', num2str(100*idx/m, '  %.01f'), ')'])
    nd=zeros(1,size(x,1));
    figure(201); clf

    subplot 312
    plot(score_vector);  xlim([1 size(x,1)]); ylim([0 1]); drawnow
    nt_mark(idx);
    ylabel('score')

    subplot 313
    colors='brgcmyk';
    hold on
    idx2=unique([idx,size(x,1)]); % add on the last sample
    old_idx=0;
    for iInterval=1:numel(idx2)
        z=nt_pca(x(old_idx+1:idx2(iInterval),:));
        dim= numel(find( mean(z.^2)/mean(z(:,1).^2) > thresh));
        nd(old_idx+1:idx2(iInterval))=dim;
        old_idx=idx2(iInterval);
    end
    z=nt_pca(x);
    nd0=numel(find( mean(z.^2)/mean(z(:,1).^2)>thresh));drawnow;
    hold on
    plot(nd0*ones(size(x,1),1), 'k--')
    plot(nd,'b'); 
    xlim([1 size(x,1)]); ylim([0 max([nd0,nd])+1]); drawnow
    nt_mark(idx);
    ylabel('ndims');
    
    subplot 311
    plot(x);  xlim([1 size(x,1)]); drawnow
    nt_mark(idx);
    if numel(idx)>1; disp(['smallest interval: ', num2str(min(diff(idx)))]); end
end

