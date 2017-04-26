function outIm = rotateImage(im,deg,newSize)
% function outIm = rotateImage(im,deg,newSize)
%
% Rotates image im by deg degrees (can be an array, in which case output is cell array)
% Click once to specify the center
% Click again to specify the direction of north wrp to the center
% Then click to specify the border of the to be rotated area wrp to the center
% (distance to center specification is taken as the radius of the new image)
% New image will be resized to size newSize*newSize
%
% By J.J.Fahrenfort, UvA 2010

% some specs
sizeIm = size(im);
middle = round(sizeIm/2);

% determine center, north and radius
figure;
imshow(im);
title('speficy center');
[ Y X ] = ginput(1);
hold on;
plot(Y,X,'ro');
%title('specify north');
%[ Y(2) X(2) ] = ginput(1);
title('specify radius');
[ Y(3) X(3) ] = ginput(1);

% compensate degrees for 'north' direction
%northAngle = atand((Y(2)-Y(1))/(X(2)-X(1)));
%deg = deg-northAngle;

% image radius to cut out
newRadius = round(sqrt((X(1)-X(3))^2 + (Y(1)-Y(3))^2));

% keep original size if no image size is specified
if nargin < 3
    newSize = 2*newRadius;
end

% make sure points fall within image
X(X<0)=0;
X(X>sizeIm(1))=sizeIm(1);
Y(Y<0)=0;
Y(Y>sizeIm(2))=sizeIm(2);

% determine how much image is needed to be able to rotate fully (add a few pixels to be sure)
pixToPad = round(sqrt(2*(newRadius^2))+10);

% pad image in all directions (mirror) so it is large enough
newIm = im;
newIm = padarray(newIm,[pixToPad pixToPad],'both','symmetric');

% determine new center location
x = X(1)+pixToPad;
y = Y(1)+pixToPad;

% cutout centerpiece
newIm = newIm(x-pixToPad:x+pixToPad,y-pixToPad:y+pixToPad,:);

% determine plot dimensions
dim1 = 3;
dim2 = 4;

% create and display new image(s)
figure;
count = 0;
for c = 1:numel(deg)
    % rotate image
    outIm{c} = imrotate(newIm,deg(c));
    % cut out new image
    newMiddle = round(size(outIm{c})/2);
    outIm{c} = outIm{c}(newMiddle(1)-newRadius:newMiddle(1)+newRadius,newMiddle(2)-newRadius:newMiddle(2)+newRadius,:);
    % resize image if necessary
    [ oldSize1 oldSize2 ~ ] = size(outIm{c});
    if all([newSize newSize] ~= [oldSize1 oldSize2])
        outIm{c} = imresize(outIm{c},[newSize newSize]);
    end
    % display image
    if mod(c,2) == 0 || c==numel(deg)
        count = count + 1;
        subplot(dim1,dim2,count);
        imshow(onlyCircle(outIm{c}));
    end
end

% If only one image, remove from cell, else make one extra plot option
if c == 1
    outIm = outIm{1};
else
    subplot(dim1,dim2,count+1);
end

% make fullscreen
set(gcf, 'Position', get(0,'Screensize'));


