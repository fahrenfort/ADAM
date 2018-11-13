clear

N=10000;

% data is concatenation of segments with different spectral structure
a=randn(N,1);
b=diff(a);
c=sin(2*pi*100*(1:N)'/N);
d=sin(2*pi*110*(1:N)'/N);
%x=[a+1;a;2*a;a;b;c];
x=[a;b;c;d];

% augment data with time-shifted versions
SHIFTS=0:10;
x=nt_multishift(x,SHIFTS);

% code 3: minimize change in covariance "shape:
figure(3); clf
DEPTH=2;    % --> 3 split points
nt_split(x,3,DEPTH);

TOL=10^-6;
disp(['rank of whole data: ', num2str(rank(x,TOL))]);        
idx=nt_split(x,3,DEPTH);
disp(['rank of intervals: ', ...
    num2str([rank(x(1:idx(1),:),TOL), ...
    rank(x(idx(1)+1:idx(2),:),TOL), ...
    rank(x(idx(2)+1:idx(3),:),TOL), ...
    rank(x(idx(3)+1:end,:),TOL)])]); 