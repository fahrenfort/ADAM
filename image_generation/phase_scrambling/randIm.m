function newIm = randIm(im)
% function newIm = randIm(im)
% shuffles all the pixels in the image and outputs the randomized image
%
% J.J.Fahrenfort, UvA 2010

imSize = size(im);
newIm = im;

for c = 1:imSize(1)
     newIm(c,:,:) = newIm(c,randperm(imSize(2)),:);
end
for c = 1:imSize(2)
     newIm(:,c,:) = newIm(randperm(imSize(1)),c,:);
end