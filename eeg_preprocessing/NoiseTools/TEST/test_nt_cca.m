
clear

x=randn(1000,11);
y=randn(1000,9);

x=nt_demean(x);
y=nt_demean(y);

[A1,B1,R1]=canoncorr(x,y);

[A2,B2,R2]=nt_cca(x,y);
A2=A2*sqrt(size(x,1));
B2=B2*sqrt(size(y,1));


figure(1); clf; 
subplot 211; plot([R1' R2']); 
if mean(A1(:).*A2(:))<0; A2=-A2; end
subplot 212; plot(([x*A1(:,1),x*A2(:,1)])); 

figure(2); clf
subplot 121; 
imagesc(nt_cov([x*A1,y*B1]));
subplot 122; 
imagesc(nt_cov([x*A2,y*B2]));

