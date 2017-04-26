function [indx1, indx2] = unpack_binned(origindx1, origindx2, oldindx1, oldindx2)
% unpacks folds on binned data into original index numbers after from
% before binning
N = size(origindx1);
for cN1 = 1:N(1)
    for cN2 = 1:N(2)
        indx1{cN1,cN2} = [oldindx1{origindx1{cN1,cN2}}]';
        indx2{cN1,cN2} = [oldindx2{origindx2{cN1,cN2}}]';
    end
end