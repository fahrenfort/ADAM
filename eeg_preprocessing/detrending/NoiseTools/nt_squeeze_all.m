function y=nt_squeeze_all(x)
%y=nt_squeeze_all(x) - squeeze structs and cell arrays
%
%  x: object to squeeze
%  
%  y: result of squeezing
%
% If x is a struct, all fields and subfields will be concatenated.  If x is
% a cell array, all elements will be concatenated.

VERBOSE=1;

x=squeeze(x);

if isnumeric (x)
    y=x;
    return;
elseif iscell(x);
    if VERBOSE; disp('cell to matrix'); end
    try
        y=cell2mat(x);
    catch
        error('could not convert cell array to matrix'); 
    end
elseif isstruct(x);
    if VERBOSE; disp('struct to cell'); end
    try
        y=struct2cell(x);
    catch
        error('could not convert structure to cell array');
    end
else
    error('!');
end

y=nt_squeeze_all(y);