clear

N=10000;

% data is concatenation of segments each of dimensionality (rank) 5
a=nt_normcol(randn(N,5)*randn(5,10));
b=nt_normcol(randn(N,5)*randn(5,10));
c=nt_normcol(randn(N,5)*randn(5,10));

% contains steps in mean, variance, covariance:
x=[a+1;a;2*a;a;b;c];

% apply nt_split with codes 1,2,3:

% code 1: minimize variance (finds steps in mean)
figure(1); clf
DEPTH=1;    % --> 1 split point
nt_split(x,1,DEPTH);

% code 2: minimize change in variance/covariance
figure(2); clf
DEPTH=2;    % --> 3 split points
nt_split(x,2,DEPTH);

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