function [M,indi,indj] = tileScramble(M,S) 
% TILESCRAMBLE - randomize tiles of a matrix
%    R = TILESCRAMBLE(M,S) randomizes the matrix M by dividing M into
%    non-overlapping blocks of the size specified by S, and shuffling these
%    blocks. M can be a N-D matrix. 
%
%    [R,I,J] = RANDBLOCK(...) also returns indices I and J, so that R equals A(I)
%    and R(J) equals A.
%
%    M can be a numerical or cell array. 
%
%  See also RAND, RANDPERM, MAT2CELL
%       and RANDSWAP, SHAKE on the File Exchange

% for Matlab R13 and up
% version 2.2 (dec 2006)
% (c) Jos van der Geest
% email: jos@jasen.nl

% History:
% Created: nov 2007
% 2.0 (dec 2007) implemented cell matrix capabilities
% 2.1 (dec 2007) added several examples for the File Exchange
% 2.2 (dec 2007) modified image example, remove minor m-lint messages

% check if the size of the blocks is specified for each dimension of M
NdimM = ndims(M);
sizM = size(M);  % size of the matrix
nS = numel(S);

% Added by J.J.Fahrenfort to simplify calling the routine for his images
if nargin < 2
    S = 8;
end
if nS < 2
    if NdimM > 2
        S = [ S S sizM(3)];
        nS = 3;
    else
        S = [ S S ];
        nS = 2;
    end
end

if nS==1,
    S = repmat(S,1,NdimM) ; % scalar expansion
elseif nS ~= NdimM,
    error('Number of elements of S should equal the number of dimensions of M.') ;
end

% S should contain positive integers only
if ~isnumeric(S) || any(S ~= fix(S)) || any(S<1),
    error('S should be a numeric array with positive integer values.') ;
end

% vectors are 2D and S can be a scalar in this case:
% the singleton dimension of M should be taken care of.
if NdimM==2 && any(sizM==1) && nS==1 && all(S~=1),
    S(sizM==1) = 1 ;    
end

if all(S==1),
    % Each element of M is a single block. We can just randomize the linear
    % indices. 
    indj = randperm(numel(M)) ;
else
    
    nb = sizM ./ S ;  % how many blocks in each dimension, this should be an integer

    B = cell(NdimM,1) ;
    % for each dimension of M: 
    for i=1:NdimM,
        % check if the number of blocks is an integer
        if nb(i) ~= fix(nb(i)), 
            error(['Size mismatch: the size of the matrix M (= %d) in dimension %d\n' ...
                    'is not an integer number of times the specified blocksize (= %d).'],...
                sizM(i),i,S(i)) ;
        else
            % expand the size in that dimension. B is used for mat2cell
            B{i} = repmat(S(i),1,nb(i)) ;
        end
    end
       
    % M can be a numerical or cell array. By randomizing the (linear)
    % indices, we do not have to worry about the actual contents of M.    
    indj = reshape(1:numel(M),sizM) ; 
    % convert to a cell array. Each cell contains a block.
    C = mat2cell(indj,B{:}) ;         
    % randomize the order of the cells
    C(randperm(numel(C))) = C ;      
    % convert back from cell array
    indj = cell2mat(C) ;              
end

% use these randomized block indices to shuffle M
M(indj) = M ;                         

if nargout>1,    
    indi(indj) = 1:numel(indj) ; 
    indi = reshape(indi,sizM) ;
    indj = reshape(indj,sizM) ;
end



