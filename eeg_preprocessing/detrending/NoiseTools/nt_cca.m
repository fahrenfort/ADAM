function [A,B,R]=nt_cca(x,y,shifts,C,m,thresh)
%[A,B,R]=nt_cca(x,y,shifts,C,m,thresh) - canonical correlation
%
%  A, B: transform matrices
%  R: r scores
%
%  x,y: column matrices
%  shifts: positive lag means y delayed relative to x
%  C: covariance matrix of [x, y]
%  m: number of columns of x
%  thresh: discard PCs below this 
%
%  Usage 1:
%   [A,B,R]=nt_cca(x,y); % CCA of x, y
%
%  Usage 2: 
%   [A,B,R]=nt_cca(x,y,shifts); % CCA of x, y for each value of shifts.
%   A positive shift indicates that y is delayed relative to x.
%
%  Usage 3:
%   C=[x,y]'*[x,y]; % covariance
%   [A,B,R]=nt_cca([],[],[],C,size(x,2)); % CCA of x,y
%
% Use the third form to handle multiple files or large data
% (covariance C can be calculated chunk-by-chunk). 
%
% C can be 3-D, which case CCA is derived independently from each page.
%
% Warning: means of x and y are NOT removed.
% Warning: A, B scaled so that (x*A)^2 and (y*B)^2 are identity matrices (differs from canoncorr).
%
% See nt_cov_lags, nt_relshift, nt_cov, nt_pca.
%
% NoiseTools.

nt_greetings; 

if ~exist('thresh','var');
    thresh=10.^-12; 
end

if exist('x','var') && ~isempty(x)
    % Calculate covariance of [x,y]
    if ~exist('y','var'); error('!'); end
    if ~exist('shifts','var')||isempty('shifts'); shifts=[0]; end
    if numel(shifts)==1 && shifts==0 && isnumeric(x) && ndims(x)==2; 
        C=[x,y]'*[x,y]; % simple case
        m=size(x,2); 
    else        
        [C,~,m]=nt_cov_lags(x,y,shifts); % lags, multiple trials, etc.
    end
    [A,B,R]=nt_cca([],[],[],C,m,thresh);
    
    if nargout==0 
        % plot something nice
        if length(shifts)>1;
            figure(1); clf;
            plot(R'); title('correlation for each CC'); xlabel('lag'); ylabel('correlation');
        end
     end
    return
end % else keep going 

if ~exist('C','var') || isempty(C) ; error('!'); end
if ~exist('m','var'); error('!'); end
if size(C,1)~=size(C,2); error('!'); end
if ~isempty(x) || ~isempty(y) || ~isempty(shifts)  ; error('!'); end
if ndims(C)>3; error('!'); end

if ndims(C) == 3
    % covariance is 3D: do a separate CCA for each page
    N=min(m,size(C,1)-m); % note that for some pages there may be fewer than N CCs
    A=zeros(m,N,size(C,3));
    B=zeros(size(C,1)-m,N,size(C,3));
    R=zeros(N,size(C,3));
    for k=1:size(C,3);
        [AA,BB,RR]=nt_cca([],[],[],C(:,:,k),m);
        A(1:size(AA,1),1:size(AA,2),k)=AA;
        B(1:size(BB,1),1:size(BB,2),k)=BB;
        R(1:size(RR,2),k)=RR;
    end
    return;
end % else keep going


%%
% Calculate CCA given C=[x,y]'*[x,y] and m=size(x,2);

% sphere x
Cx=C(1:m,1:m);
[V, S] = eig(Cx) ;  
V=real(V); S=real(S);
[E, idx] = sort(diag(S)', 'descend') ;
keep=find(E/max(E)>thresh);
topcs = V(:,idx(keep));
E = E (keep);
EXP=1-10^-12; 
E=E.^EXP; % break symmetry when x and y perfectly correlated (otherwise cols of x*A and y*B are not orthogonal)
A1=topcs*diag(sqrt((1./E)));

% sphere y
Cy=C(m+1:end,m+1:end);
[V, S] = eig(Cy) ;  
V=real(V); S=real(S);
[E, idx] = sort(diag(S)', 'descend') ;
keep=find(E/max(E)>thresh);
topcs = V(:,idx(keep));
E = E (keep);
E=E.^EXP; % 
A2=topcs*diag(sqrt((1./E)));

% apply sphering matrices to C
AA=zeros( size(A1,1)+size(A2,1), size(A1,2)+size(A2,2) );
AA( 1:size(A1,1), 1:size(A1,2) )=A1;
AA( size(A1,1)+1:end, size(A1,2)+1:end )=A2;
C= AA' * C * AA;

N=min(size(A1,2),size(A2,2)); % number of canonical components

% PCA
[V, S] = eig(C) ;
%[V, S] = eigs(C,N) ; % not faster
V=real(V); S=real(S);
[E, idx] = sort(diag(S)', 'descend') ;
topcs = V(:,idx);

A=A1*topcs(1:size(A1,2),1:N)*sqrt(2);  % why sqrt(2)?...
B=A2*topcs(size(A1,2)+1:end,1:N)*sqrt(2);
R=E(1:N)-1; 


%{
Why does it work?
If x and y were uncorrelated, eigenvalues E would be all ones. 
Correlated dimensions (the canonical correlates) should give values E>1, 
i.e. they should map to the first PCs. 
To obtain CCs we just select the first N PCs. 
%}

%%

%%
% test code
if 0
    % basic
    clear
    x=randn(10000,20);
    y=randn(10000,8);
    y(:,1:2)=x(:,1:2); % perfectly correlated
    y(:,3:4)=x(:,3:4)+randn(10000,2); % 1/2 correlated
    y(:,5:6)=x(:,5:6)+randn(10000,2)*3; % 1/4 correlated
    y(:,7:8)=randn(10000,2); % uncorrelated
    [A,B,R]=nt_cca(x,y);
    figure(1); clf
    subplot 321; imagesc(A); title('A');
    subplot 322; imagesc(B); title('B');
    subplot 323; plot(R, '.-'); title('R')
    subplot 324; nt_imagescc((x*A)'*(x*A)); title ('covariance of x*A');
    subplot 325; nt_imagescc((y*B)'*(y*B)); title ('covariance of y*B');
    subplot 326; nt_imagescc([x*A,y*B]'*[x*A,y*B]); title ('covariance of [x*A,y*B]');
end

if 0 
    % compare with canoncorr
    clear
    x=randn(1000,11);
    y=randn(1000,9);
    x=x-repmat(mean(x),size(x,1),1); % center, otherwise result may differ slightly from canoncorr
    y=y-repmat(mean(y),size(y,1),1);
    [A1,B1,R1]=canoncorr(x,y);
    [A2,B2,R2]=nt_cca(x,y);   
    A2=A2*sqrt(size(x,1)); % scale like canoncorr
    B2=B2*sqrt(size(y,1));
    figure(1); clf; 
    subplot 211; 
    plot([R1' R2']); title('R'); legend({'canoncorr', 'nt_cca'}, 'Interpreter','none'); 
    if mean(A1(:,1).*A2(:,1))<0; A2=-A2; end
    subplot 212; 
    plot(([x*A1(:,1),x*A2(:,1)])); title('first component'); legend({'canoncorr', 'nt_cca'}, 'Interpreter','none'); 
    figure(2); clf;set(gcf,'defaulttextinterpreter','none')
    subplot 121; 
    nt_imagescc([x*A1,y*B1]'*[x*A1,y*B1]); title('canoncorr'); 
    subplot 122; 
    nt_imagescc([x*A2,y*B2]'*[x*A2,y*B2]); title('nt_cca');
end

if 0
    % time
    x=randn(100000,100); 
    tic; 
    [A,B,R]=nt_cca(x,x); 
    disp('nt_cca time: ');
    toc    
    [A,B,R]=canoncorr(x,x); 
    disp('canoncorr time: ');
    toc
%     [A,B,R]=cca(x,x); 
%     disp('cca time: ');
%     toc
end

if 0
    % lags
    x=randn(1000,10);
    y=randn(1000,10);
    y(:,1:3)=x(:,1:3);
    lags=-10:10;
    [A1,B1,R1]=nt_cca(x,y,lags);
    figure(1); clf
    plot(lags,R1'); xlabel('lag'); ylabel('R');
end

if 0
    % what happens if x & y perfectly correlated?
    x=randn(1000,10);
    y=randn(1000,10); y=x(:,randperm(10)); %+0.000001*y;
    [A1,B1,R1]=nt_cca(x,y);
    figure(1); clf
    nt_imagescc([x*A1,y*B1]'*[x*A1,y*B1]);
end    

if 0
    % x and y are cell arrays
    x=randn(1000,10); 
    y=randn(1000,10);
    xx={x,x,x};  yy={x,y,y};
    [A,B,R]=nt_cca(xx,yy);
    disp('seems to work...');
end

    

