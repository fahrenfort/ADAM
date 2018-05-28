function capfile = trycapfile
% try to find the capfile
try 
    capfile = findcapfile;
    if ~exist(capfile,'file')
        error('uggh');
    end
catch
    disp('Cannot find the hardcoded capfile location, trying dynamically');
    capfile = which('standard-10-5-cap385.elp');
end