function [C,A,score]=nt_cluster1D_b(x);
%[C,A,score]=nt_cluster1D_b(x) - cluster 1D data into 2 clusters
%
%  x: column vector or matrix of data to cluster
%
%  C: centroid pairs (one pair per column)
%  A: ownership matrix (0, 1)
%  score: energy/total energy, for each column

if nargin<1; error('!'); end
if size(x,1)<2; error('too small to cluster'); end

A=zeros(size(x));       % cluster ownership labels
C=zeros(2,size(x,2));   % centroids

for iCol=1:size(x,2)
    
    xx=x(:,iCol);
    [xx,iSort]=sort(xx);
    [idx,score_vector,score0]=nt_split(xx);
    score(:,iCol)=score0;
    C(:,iCol)=[mean(xx(1:idx)),mean(xx(idx+1:end))];
    t=1:size(xx,1);
    A(t(iSort(idx+1:end)), iCol)=1; % 0: low values, 1: high values
    
%     figure(1); clf; subplot 211;
%     hold on; histogram(xx,-5:0.01:5, 'displaystyle','stairs' ); histogram(xx(1:idx),-5:0.01:5); plot(C(:,iCol),[500 500], '.r');
%     subplot 212; plot(x(:,iCol))
%     disp(score0)
%     pause
end


