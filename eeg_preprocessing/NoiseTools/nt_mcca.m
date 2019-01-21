function [A,score,AA]=nt_mcca(C,N,Nkeep)
%[A,score,AA]=nt_mcca(C,N) - multiple cca
%
%  A: transform matrix
%  score: commonality score (ranges from 1 to N)
%  AA: array of subject-specific MCCA transform matrices
% 
%  C: covariance matrix of aggregated data sets
%  N: number of channels of each data set
%  Nkeep: number of components to keep (for orthogonal transforms)

if nargin<3; Nkeep=[]; end
if nargin<2; error('!'); end
if size(C,1) ~= size(C,2); error('!'); end
if size(C,1) ~= round(size(C,1)/N)*N; error('!'); end

% sphere by blocks
nblocks=size(C,1)/N;
for iBlock=1:nblocks
    idx=(iBlock-1)*N + (1:N);
    CC=C(idx,idx);
    [V, S] = eig(CC) ;  
    V=real(V); S=real(S);
    [E,idx2] = sort(diag(S)', 'descend');
    topcs=V(:,idx2);
    EXP=1-10^-12; 
    E=E.^EXP; % break symmetry when x and y perfectly correlated (otherwise cols of x*A and y*B are not orthogonal)
    EE=(1./E); EE(find(E<=0))=0;
    A(idx,idx)=topcs*diag(sqrt(EE));
end
C=A'*C*A;


% final PCA
[V, S] = eig(C) ;
V=real(V); S=real(S);
[E, idx] = sort(diag(S)', 'descend') ;
topcs = V(:,idx);
A=A*topcs;
%A=A(:,1:N);

C=topcs'*C*topcs;
score=diag(C);


if nargout>2;
    AA=[];
    for iBlock=1:nblocks
        AA{iBlock}=A(N*(iBlock-1)+(1:N),:);
    end
end

if nargout>3;
    if isempty(Nkeep); error('must specify Nkeep'); end
    AAA=[];
    for iBlock=1:nblocks
        % covariance of subject's data
        idx=(iBlock-1)*N + (1:N);
        C11=C(idx,idx);
        % covariance of selected MCCA components
        tmp=A(:,1:Nkeep);
        C22=tmp'*C*tmp;
        % cross covariance between subject's data and transformed data
        C12=C(idx,:)*tmp; clear tmp
        C21=C12';
        % CCA:
        [tmp]=nt_cca([],[],[],[C11,C12;C21,C22],N);
        AAA{iBlock}=tmp;
    end
end


return

% test code

% 3 uncorrelated data sets
figure(1); clf;
x1=randn(10000,10); x2=randn(10000,10); x3=randn(10000,10); 
x=[x1,x2,x3];
C=x'*x;
[A,score,AA]=nt_mcca(C,10);
subplot 131; nt_imagescc(A); title('mcca transform');
subplot 132; nt_imagescc(A'*C*A); title('covariance of transformed data');
subplot 133; nt_imagescc(x'*(x*A)); title('crosscorr between raw & transf'); xlabel('transformed'); ylabel('raw');
z=x*A;
figure(11); clf;
plot(mean(z.^2));

% 3 identical data sets
figure(2); clf
x1=randn(10000,10); x=[x1,x1,x1]; 
C=x'*x; 
[A,score,AA]=nt_mcca(C,10);
subplot 131; nt_imagescc(A); title('mcca transform');
subplot 132; nt_imagescc(A'*C*A); title('covariance of transformed data');
subplot 133; nt_imagescc(x'*(x*A)); title('cross correlation'); xlabel('transformed'); ylabel('raw');

% 3 data sets with shared parts
figure(3); clf
x1=randn(10000,5); x2=randn(10000,5); x3=randn(10000,5); x4=randn(10000,5); 
x=[x2,x1,x3,x1,x4,x1];
C=x'*x; 
[A,score,AA]=nt_mcca(C,10);
subplot 131; nt_imagescc(A); title('mcca transform');
subplot 132; nt_imagescc(A'*C*A); title('covariance of transformed data');
subplot 133; nt_imagescc(x'*(x*A)); title('cross correlation'); xlabel('transformed'); ylabel('raw');

% 3 data sets with parts shared 2 by 2
figure(4); clf
x1=randn(10000,5); x2=randn(10000,5); x3=randn(10000,5); x4=randn(10000,5); 
x=[x2,x1,x3,x1,x3,x4];
C=x'*x; 
[A,score,AA]=nt_mcca(C,10);
subplot 131; nt_imagescc(A); title('mcca transform');
subplot 132; nt_imagescc(A'*C*A); title('covariance of transformed data');
subplot 133; nt_imagescc(x'*(x*A)); title('cross correlation'); xlabel('transformed'); ylabel('raw');

