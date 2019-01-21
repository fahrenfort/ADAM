function r=nt_sns0(c,nneighbors,skip,wc)
% r=nt_sns0(c,nneigbors,skip,wc) - sensor noise suppression
%
%   r: denoising matrix
%
%   c: full covariance of data to denoise
%   nneighbors: number of channels to use in projection 
%   skip: number of neighbors to skip [default: 0]
%   wc: weighted covariance
%
% 


n=size(c,1);

if nargin<2 || isempty(nneighbors); error('need to specify nneighbors'); end
if nargin<3 || isempty(skip); skip=0; end
if nargin<4 || isempty(wc); wc=c; end

nneighbors=min(nneighbors,n-skip-1);

r=zeros(size(c));

% normalize
d=sqrt(1./(diag(c)+eps));
c=nt_vecmult(nt_vecmult(c,d),d');

for iChan=1:n
 
    c1=c(:,iChan);                      % correlation of channel with all other channels
    [c1,idx]=sort(c1.^2,1,'descend');   % sort by correlation
    idx=idx(skip+2:skip+1+nneighbors);  % keep best

    % pca neighbors to orthogonalize them
    c2=wc(idx,idx);
    [topcs,eigenvalues]=nt_pcarot(c2);
    topcs=topcs*diag(1./sqrt(eigenvalues));
    topcs(find(isinf(topcs)|isnan(topcs)))=0;
    
    % augment rotation matrix to include this channel
    topcs=[1,zeros(1,nneighbors);zeros(nneighbors,1),topcs];
    
    % correlation matrix for rotated data
    c3=topcs'*wc([iChan;idx],[iChan;idx])*topcs;
    
    % first row defines projection to clean component iChan
    c4=c3(1,2:end)*topcs(2:end,2:end)';

    % insert new column into denoising matrix
    r(idx,iChan)=c4;

end
