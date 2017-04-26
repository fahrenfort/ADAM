function measure = tuning_from_matrix(ctfMatrix,method,crossclass)
% compute single tuning measure from a ctfMatrix
% J.J.Fahrenfort, VU, 2016

if nargin < 3
    crossclass = true;
end
if nargin < 2
    method = 'slope';
end

if isvector(ctfMatrix) % does not contain time
    measure = tuning_from_matrix_only(ctfMatrix,method);
elseif ndims(ctfMatrix) == 3 % wrapper in case ctfMatrix also contains time as the first two dimensions
    measure = zeros(size(ctfMatrix,1), size(ctfMatrix,2));
    t2start = 1;
    t2stop = size(ctfMatrix,2);
    for t=1:size(ctfMatrix,1) % test loop
        if ~crossclass % train loop (or diagonal)
            t2start = t;
            t2stop = t;
        end
        for t2=t2start:t2stop
            measure(t,t2) = tuning_from_matrix_only(squeeze(ctfMatrix(t,t2,:)),method);
        end
    end
else
    error('the ctfMatrix array does not contain a sensible number of dimensions');    
end

function measure = tuning_from_matrix_only(ctfMatrix,method)
% compute accuracy measure
if strcmpi(method,'slope')
    [ ~, ~, measure ] = fit_slope(ctfMatrix');
elseif sum(strcmpi(method,{'ctf_sigma','ctf_A'}))
    [ ~, ~, sigma, A ] = fit_ctf(ctfMatrix');
    if strcmpi(method,'ctf_sigma')
        measure = sigma;
    elseif strcmpi(method,'ctf_A')
        measure = A;
    end
else
    error('no existing method was specified');
end