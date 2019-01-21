function [z,idx]=nt_tsregress(x,y,shifts,xw,yw,keep,threshold)
%[z,idx]=nt_tsregress(x,y,shifts,xw,yw,keep,threshold) - time-shift regression
%
%  z: part of x modeled by time-shifted y
%  idx: x(idx) maps to z
%
%  x: data to model
%  y: regressor
%  shifts: array of shifts to apply (default: [0])
%  xw: weights to apply to x
%  yw: weights to apply to y
%  keep: number of components of shifted regressor PCs to keep (default: all)
%  threshold: discard PCs with eigenvalues below this (default: 0)
%
% Data X are regressed on time-shifted versions of Y. X and Y are initially 
% time-aligned, but because of the shifts, Z is shorter than X.  Z is
% time-aligned with X(IDX).

if nargin<2; error('!'); end
if nargin<3||isempty(shifts); shifts=[0]; end
if nargin<4; xw=[]; end
if nargin<5; yw=[]; end
if nargin<6; keep=[]; end
if nargin<7; threshold=[]; end

if size(x,1) ~= size(y,1); error('!'); end

% shifts must be non-negative
mn=min(shifts);
if mn<0; 
    shifts=shifts-mn; 
    x=x(-mn+1:end,:,:);
    y=y(-mn+1:end,:,:);
end
nshifts=numel(shifts);

% % flag outliers in x and y
% if ~isempty(toobig1) || ~isempty(toobig2)
%     xw=nt_find_outliers(x,toobig1,toobig2);
%     yw=nt_find_outliers(y,toobig1,toobig2);
% else
%     xw=[];yw=[];
%     %xw=ones(size(x)); yw=ones(size(y));
% end

% subtract weighted means

if ndims(x)==3    
    [Mx,Nx,Ox]=size(x);
    [My,Ny,Oy]=size(y);
    x=nt_unfold(x);
    y=nt_unfold(y);
    [x,xmn]=nt_demean(x,xw);
    [y,ymn]=nt_demean(y,yw);
    x=nt_fold(x,Mx);
    y=nt_fold(y,My);
else
    [x,xmn]=nt_demean(x,xw);
    [y,ymn]=nt_demean(y,yw);
end


% covariance of y
[cyy,totalweight]=nt_cov(y,shifts',yw);
cyy=cyy./totalweight;

% cross-covariance of x and y
[cxy, totalweight]=nt_cov2(x,y,shifts',xw,yw);
disp('!!!!!!!!!   WARNING: calling obsolete code  !!!!!!!!!!!!!!!!');
%[cxy, totalweight]=nt_xcov(x,y,shifts',xw,yw);
cxy=cxy./totalweight;

% regression matrix
r=nt_regcov(cxy,cyy,keep,threshold);
    
% regression
if ndims(x)==3
    x=nt_unfold(x);
    y=nt_unfold(y);
    [m,n]=size(x);
    mm=m-max(shifts);
    z=zeros(size(x));
    for k=1:nshifts
        kk=shifts(k);
        idx1=kk+1:kk+mm;
        idx2=k+(0:size(y,2)-1)*nshifts;
        z(1:mm,:)=z(1:mm,:)+y(idx1,:)*r(idx2,:);
    end
    z=nt_fold(z,Mx);
    z=z(1:end-max(shifts),:,:);
else
    [m,n]=size(x);
    z=zeros(m-max(shifts),n);
    for k=1:nshifts
        kk=shifts(k);
        idx1=kk+1:kk+size(z,1);
        %idx2=k*size(y,2)+1:(k+1)*size(y,2);
        idx2=k+(0:size(y,2)-1)*nshifts;
        z=z+y(idx1,:)*r(idx2,:);
    end
end

% idx allows x to be aligned with z
offset=max(0,-mn);
idx=offset+1:offset+size(z,1);

%elseif ndims(x)==3
    
%end

