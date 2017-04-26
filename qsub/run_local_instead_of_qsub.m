function run_local_instead_of_qsub(path_on_lisa,function_name, settings, varargin)
% just to check, takes the same arguments but runs locally
% local_path_of_function contains the path where the function is on locally e.g.
% '/home/johannes/WM_EEG/mvpa_scripts'
% function_name contains the function name, e.g.
% 'classify_eeg'
% settings.distantpath will be replaced with
% settings.localpath 
% varargin is a series of arguments the function may take
% each element can be a cell array, for which all combination will run
% walltime is disregarded
% J.J.Fahrenfort, VU, 2014, 2015
% we need to create a loop for each var in varargin
combMat = allcombs(varargin{:});

% run each combination of arguments
for cQsubs = 1:size(combMat,1)
    if exist(function_name,'file')
        line = sprintf('%s(',function_name);
    else
        line = sprintf('%s%s%s(',path_on_lisa,filesep,function_name);
    end
    for cArgs = 1:size(combMat,2)
        line = [line sprintf('''%s'',',combMat{cQsubs,cArgs})];
    end
    line = [line(1:end-1) ');'];
    % replace with local directory
    if ~isempty(settings)
        line = strrep(line,settings.distantpath,settings.localpath);
    end
    % execute
    disp(line);
    eval(line);
end
return

function combMat = allcombs(varargin)
for c = 1:nargin
    if isnumeric(varargin{c})
        varargin{c} = strsplit(num2str(varargin{c}),' ');
    end
    if ischar(varargin{c})
        varargin{c} = {varargin{c}};
    end
end
sizeVec = cellfun('prodofsize', varargin);
indices = fliplr(arrayfun(@(n) {1:n}, sizeVec));
[indices{:}] = ndgrid(indices{:});
combMat = cellfun(@(c,i) {reshape(c(i(:)), [], 1)}, ...
    varargin, fliplr(indices));
combMat = [combMat{:}];
return