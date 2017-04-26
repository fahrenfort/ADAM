function outIm = onlyCircle(im, radius)
% function outIm = onlyCircle(im, radius)
% turns the part of the image im that is not in a circle with radius radius
% white and returns the resulting image.
%
% By J.J.Fahrenfort, UvA 2010

imSize = size(im);

if nargin <2
     radius = min(imSize(1:2))/2;
end

% reshape to rgb array and normalize
if numel(imSize)<3
    % double(max(im(:)))
    outIm = repmat(double(im)./double(max(im(:))),[1 1 3]);
end

% change all that is outside circle to green
for x = 1:imSize(1)
    for y = 1:imSize(2)
        if sqrt((x-imSize(1)/2)^2 + (y-imSize(2)/2)^2) > radius
            outIm(x,y,:) = [0 255 0] ;
        end
    end
end

