function [c,tw]=nt_cov(x,shifts,w);
%[c,tw]=nt_cov(x,shifts,w) - time shift covariance
%
%  c: covariance matrix
%  tw: total weight (c/tw is normalized covariance)
%
%  x: data
%  shifts: array of time shifts (must be non-negative)
%  w: weights
%  
% If SHIFTS==[0], this function calculates the covariance matrix of X, 
% weighted by W.  In the general case it calculates the covariance matrix
% of time-shifted X.
%
% X can be 1D, 2D or 3D.  It can also be a cell array. 
% W can be 1D (if X is 1D or 2D) or 2D (if X is 3D). The same weight is
% applied to each column.
% 
% Output is a 2D matrix with dimensions (ncols(X)*numel(SHIFTS))^2.
% It is made up of an ncols(X)*ncols(X) matrix of submatrices, each of 
% dimensions numel(SHIFTS)*numel(SHIFTS).
%
% NoiseTools

%% arguments
if nargin<3; w=[]; end;
if nargin<2||isempty(shifts); shifts=0; end;
if prod(size(x))==0; error('data empty'); end
shifts=shifts(:);           % --> column vector
nshifts=numel(shifts); 


if isempty(w) 
    %% no weights
    if isnumeric(x)
        %% matrix
        [m,n,o]=size(x);
        c=zeros(n*nshifts);
        for k=1:o
            xx=nt_multishift(x(:,:,k),shifts);
            c=c+xx'*xx;
        end
        tw=size(xx,1)*o;
    elseif iscell(x)
        %% cell array
        [m,n]=size(x{1});
        o=length(x);
        c=zeros(n*nshifts);
        for k=1:o;
            if size(x{k},2)~=n;  disp([size(x{k}),n, k]); error('!'); end
            xx=nt_multishift(x{k}(:,:),shifts);
            c=c+xx'*xx;
        end
        tw=size(xx,1)*o;
    else error('!'); end   
else
    %% weights
    if isnumeric(x)
        %% matrix
        [m,n,o]=size(x);
        c=zeros(n*nshifts);
        for k=1:o
            if ~all(shifts == [0]); 
                xx=nt_multishift(x(:,:,k),shifts); 
                ww=nt_multishift(w(:,:,k),shifts);
                ww=min(ww,[],2);
            else
                xx=x(:,:,k); ww=w(:,:,k);
            end
            xx=nt_vecmult(xx,ww);
            c=c+xx'*xx;
        end
        tw=sum(w(:));
    else
        %% cell array
        c=zeros(n*nshifts);
        [m,n]=size(x{1});
        o=length(x);
        for k=1:o
            if ~all(shifts == [0]); 
                xx=nt_multishift(x{k}(:,:),shifts); 
                ww=nt_multishift(w{k}(:,:),shifts);
                ww=min(ww,[],2);
            else
                xx=x{k}(:,:); ww=w{k}(:,:);
            end
            xx=nt_vecmult(xx,ww);
            c=c+xx'*xx;
        end
        tw=sum(w(:));
    end       
end

