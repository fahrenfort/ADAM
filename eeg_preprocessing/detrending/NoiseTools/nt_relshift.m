function [xx,yy]=nt_relshift(x,y,shift)
%[xx,yy]=nt_relshift(x,y,shift,flag) - delay x relative to y 
%
%  xx, yy: shifted matrices
%
%  x,y: column matrices to shift
%  shift: amount to delay x (can be negative or fractionary)
%  
% If shift has multiple values, xx and yy are 3D matrices, one shift per
% page. Shifted data are zero-padded.
%
% NoiseTools

if nargin<3; error('!'); end
if ~isnumeric(x); error('!'); end
if size(x,1)~=size(y,1); 
%    warning(['x and y have different nrows: ', num2str([size(x,1), size(y,1)])]);
    m=min(size(x,1),size(y,1));
    x=x(1:m,:,:); 
    y=y(1:m,:,:);
    %error('!'); 
end

if shift ~= round(shift); error('fractionary shifts not yet implemented'); end

if length(shift)==1
    if shift>0
        yy=y(1:end-shift,:);
        xx=x(shift+1:end,:);
    else
        yy=y(-shift+1:end,:);
        xx=x(1:end+shift,:);
    end   
else
    xx=zeros(size(x,1), size(x,2), length(shift));    
    yy=zeros(size(y,1), size(y,2), length(shift));
    for iShift=1:length(shift)
        s=shift(iShift);
        if s>0
            yy(1:end-s,:,iShift)=y(1:end-s,:);
            xx(1:end-s,:,iShift)=x(s+1:end,:);
        else
            yy(1:end+s,:,iShift)=y(-s+1:end,:);
            xx(1:end+s,:,iShift)=x(1:end+s,:);
        end   
    end
end

if 0 
    x=sin(2*pi*3*(1:1000)'/1000);
    y=x;
    figure(1); clf;
    subplot 131;
    [xx,yy]=nt_relshift(x,y,100);
    plot([xx,yy])
    subplot 132; 
    [xx,yy]=nt_relshift(x,y,-100:10:100);
    plot(squeeze(xx));
    subplot 133; 
    plot(squeeze(yy));
end
    
    
