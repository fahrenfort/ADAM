function nt_imagescc(C)
%nt_imagescc - plot image with symmetric scaling

m=max(abs(C(:)));
imagesc(C,[-m-realmin,m+realmin]);
