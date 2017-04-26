function [logstruct] = readPsychoPy(logdatadir,logfiles,sessions)
% function [logstruct] = readPsychoPy(logdatadir,logfiles,sessions)
if nargin<2
    logfiles = dir([logdatadir filesep '*.csv']);
    logfiles = {logfiles.name};
end
longestCommon = LCS(logfiles{:});
for cSes = 1:numel(logfiles)
    % logfile
    [~,logfile,~] = fileparts(logfiles{cSes});
    % session name, keep it as short as possible
    if nargin<3
        session = strrep(logfile,longestCommon,'_');
        if ~isnan(str2double(logfile(1)))
            session = ['S' session];
        end
    else
        session = sessions{cSes};
    end
    logfile = [logdatadir filesep logfile '.csv'];
    % read data
    data = readtable(logfile);
    % little hack to get the session and task from the table
    if any(strcmp('session',fieldnames(data))) && any(strcmp('task_type',fieldnames(data)))
        if numel(unique(data.task_type)) == 1 && numel(unique(data.session)) == 1
            session = [ session(1:min(strfind(session,'_')-1)) '_S0' num2str(unique(data.session)) '_B' num2str(unique(data.task_type)) ];
            %session = [ session(1:min(strfind(session,'_')-1)) '_ses' num2str(unique(data.session)) '_task' num2str(unique(data.task_type)) ];
        end
    end
    logstruct.(session) = data;
end
logstruct = orderfields(logstruct);

