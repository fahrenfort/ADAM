% Apply MCCA to Sam's fMRI data
%
% MCCA boils down to: (1) PCA, normalize each subject, (2) concatenate
% along the PC dimension, (3) PCA.

clear;

path_to_data_file = '/data/meg/sam/for-alain-3mm.mat';
load(path_to_data_file, 'X');
load(path_to_data_file, 'category_labels');
[~,iSortCategories]=sort(category_labels);

% right hemisphere upside-down?
X{1}=flipud(X{1});
% flip both hemispheres left-right
X{1}=fliplr(X{1});
X{2}=fliplr(X{2});

if 1
    % scale data to remove means and normalize
    flags=[1 1 0 0]; % for each pixel remove mean and normalize
    X=scale_data(X,flags);
end

[~,~,nSnd,nRep,nSubj]=size(X{1});

activations_all=[];
for iSubj=1:nSubj  
    
    % merge data for both hemispheres
    xx=[];
    for iHem=1:2
        x=X{iHem}(:,:,:,:,iSubj);
        x(find(isnan(x)))=0;
        [nX,nY,~,~,~]=size(x);
        x=reshape(x,nX*nY,nSnd,nRep);
        xx=[xx;x];      % concatenate unfolded pixels
    end
    xx=mean(xx,3);    % average over repeats
    xx(find(isnan(xx)))=0;
    % --> xx has size npixels*nsounds
    
    xx=xx(:,randperm(165));
    
    % more pixels than sounds --> PCA over the pixel dimension
    xx=nt_demean(xx);
    topcs=nt_pca0(xx);   % each PC defined by an activation pattern over sounds, and a pixel pattern
    frompcs=pinv(topcs); % columns of frompcs are activation patterns
    NPCs=15;
    activations_all=[activations_all,frompcs(1:NPCs,:)']; % concatenate    
end

% mCCA: normalize concatenated PCs then PCA
activations_all=nt_normcol(activations_all);
z=nt_pca(activations_all);

figure(1); clf; set(gcf,'position', [440 378 600 200], 'color',[1 1 1]);

axes('position', [.32 .2 .16 .70])
plot((mean(z.^2)),'.-k','markersize', 22); 
xlim([0 size(z,2)]); ylim([0 11])
set(gca,'ytick',[0 5 10], 'fontsize',14); xlabel('SC'); %ylabel('variance');
line([0 size(z,1)],[1 1],'linestyle',':');  title('(b)');


activations_all=[];
for iSubj=1:nSubj  
    
    % merge data for both hemispheres
    xx=[];
    for iHem=1:2
        x=X{iHem}(:,:,:,:,iSubj);
        x(find(isnan(x)))=0;
        [nX,nY,~,~,~]=size(x);
        x=reshape(x,nX*nY,nSnd,nRep);
        xx=[xx;x];      % concatenate unfolded pixels
    end
    xx=mean(xx,3);    % average over repeats
    xx(find(isnan(xx)))=0;
    % --> xx has size npixels*nsounds
    
    % more pixels than sounds --> PCA over the pixel dimension
    xx=nt_demean(xx);
    topcs=nt_pca0(xx);   % each PC defined by an activation pattern over sounds, and a pixel pattern
    frompcs=pinv(topcs); % columns of frompcs are activation patterns
    NPCs=15;
    activations_all=[activations_all,frompcs(1:NPCs,:)']; % concatenate    
end

% mCCA: normalize concatenated PCs then PCA
activations_all=nt_normcol(activations_all);
z=nt_pca(activations_all);

axes('position', [.54 .2 .16 .70])
plot((mean(z.^2)),'.-k','markersize', 22); 
xlim([0 size(z,2)]); ylim([0 11])
set(gca,'ytick',[0 5 10], 'fontsize',14); xlabel('SC'); %ylabel('variance');
line([0 size(z,1)],[1 1],'linestyle',':');  title('(c)');

xx=[];
for iSet=1:10
    xx=[xx, nt_normcol(nt_pca(randn(100000,15)))];
end
z=nt_pca(xx);

axes('position', [.1 .2 .16 .70])
plot((mean(z.^2)),'.-k','markersize', 22); 
xlim([0 size(z,2)]); ylim([0 11])
set(gca,'ytick',[0 5 10], 'fontsize',14); xlabel('SC'); ylabel('variance');
line([0 size(z,1)],[1 1],'linestyle',':');  title('(a)');

      
%subplot 144;
axes('position', [.76 .2 .16 .70])
load Fig_A_data
plot(p,'.-k','markersize', 22); 
xlim([0 numel(p)]); ylim([0 11])
set(gca,'ytick',[0 5 10], 'fontsize',14); xlabel('SC'); %ylabel('variance');
line([0 size(z,1)],[1 1],'linestyle',':'); title('(d)');




set(gcf, 'PaperPositionMode', 'auto');
print ('-depsc2', 'PAPER_MCCA4/Fig_variance_profiles')


function [Y,a]=scale_data(X,flags)
%[Y,a]=scale_data(X,flags) - rescale data over pixels
% 
%  Y: scaled data
%  a: scaling data
%
%  X: data to scale
%  flags

if nargin<2 || isempty(flags); flags=[0 0 0 0]; end
if numel(flags<4); flags(4)=0; end


for iHemi=1:2
    [nX,nY,nSnd,nRep,nSubj] = size(X{iHemi});
    for iSubj=1:nSubj;
        x=X{iHemi}(:,:,:,:,iSubj);
        iNan=find(isnan(x));
        x(iNan)=0;
        y=x;

        if flags(1) % remove mean for each pixel
            m=mean(x(:,:,:),3);
            y=bsxfun(@minus,y,m);
            a{1,iSubj}=m;
        end

        if flags(2) % normalize each pixel
            n=sqrt(mean(y(:,:,:).^2,3));
            y=bsxfun(@times,y,1./n);
            y(find(isnan(y)))=0;
            a{2,iSubj}=sqrt(n);
        end

        if flags(3) % remove mean for each repeat
            m=mean(reshape(x,nX*nY*nSnd,nRep),1);
            m=reshape(m,[1 1 1 numel(m)]);
            x=bsxfun(@minus, x, m);
            a{3,iSubj}=m;
        end
        
        if flags(4) % normalize each repeat
            n=sqrt(mean(reshape(x.^2,nX*nY*nSnd,nRep),1));
            n=reshape(n,[1 1 1 numel(n)]);
            [~,~,nSnd,nRep,nSubj] = size(X{1});
            y(find(isnan(y)))=0;
            a{4,iSubj}=sqrt(n);
        end
            
        y(iNan)=nan;
        Y{iHemi}(:,:,:,:,iSubj)=y;
    end
end
end
