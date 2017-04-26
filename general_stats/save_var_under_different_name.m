function save_var_under_different_name(fullfilename,varargin)
% function save_var_under_different_name(fullfilename,varargin)
% Saves variables in a workspace under a different name in file
% fullfilename. Only works on existing files, always appends. Waits for
% file creation before starting.
% use fullfilename for file name and pairs of variables containing data and 'newvarname' 
% example:
% save_var_under_different_name('dirname\my_matfile.mat',data1,'dynamicdataname1',data2,'dynamicdataname2')
% only work when appending
% J.J.Fahrenfort, VU 2016
if mod(numel(varargin),2)
    error('varargin must be even for save_var_under_different_name to work');
end
for c = 1:2:numel(varargin)
    if ~isempty(varargin{c})
        S.(varargin{c+1}) = varargin{c};
    end
end
% make sure it ends in .mat
fullfilename = regexprep(fullfilename,'.mat','');
fullfilename = [fullfilename '.mat'];
c = 0;
while ~exist(fullfilename,'file')
    c = c + 1;
    disp(['waiting for ' fullfilename ' to be created.']);
    pause(1);
    if c == 120
        error([ fullfilename ' does not exist, so cannot append.']);
    end
end
% it exists, so append
save(fullfilename, '-struct', 'S', '-v7.3', '-append');
