%% scripts to read trigger values
clear
clc
warning('off','all');

eegdatadir = '/Users/VU-MBP/Dropbox/Work/Experimenten/WM/BDFs';
logdatadir = '/Users/VU-MBP/Dropbox/Work/Experimenten/WM/logfiles';

%
% '01_WM_2014_jul_09_0958' 
sessions = { 'WM_pilot_pp1s1' 'WM_pilot_pp1s2' 'WM_pilot_pp1s2_2' 'WM_pilot_pp2s1' 'WM_pilot_pp2s2' 'WM_pilot_pp3s1' 'WM_pilot_pp3s2' 'WM_pilot_pp4s1' 'WM_pilot_pp4s2'   } ;
logfiles = { '01_WM_2014_jul_08_1009' '01_WM_2014_jul_09_0958' '01_WM_2014_jul_09_0958' '02_WM_2014_jul_08_1229' '02_WM_2014_jul_09_1207' '03_WM_2014_jul_08_1521' '03_WM_2014_jul_09_1458' '04_WM_2014_jul_15_1004' '04_WM_2014_jul_16_1031' };

for cSes = 1:numel(sessions)
    % session name
    session = sessions{cSes};
    log = logfiles{cSes};
    
    % get filenames
    eegfile = [eegdatadir filesep session '.bdf'];
    logfile = [logdatadir filesep log '.csv'];
    % read data
    [ triggers.(session), logtable.(session) ] = jjf_triggers_from_logs(eegfile,logfile);
    behavior(cSes) = jjf_computeBehavior(logtable.(session),session);
    disp(behavior(cSes))
end
save([eegdatadir filesep 'triggercodes_new'], 'triggers');
%save([eegdatadir filesep 'behavior'], 'behavior', 'logtable');
warning('on','all');

%% plot subjects 
subjects = fields(logtable);
nSubj = numel(subjects);
barH = figure;
retH = figure;
for cSubj = 1:nSubj;
    figure(barH);
    subplot(sqrt(nSubj),sqrt(nSubj),cSubj);
    hist(logtable.(subjects{cSubj}).bar_duration(logtable.(subjects{cSubj}).setsize==10));
    title(strrep(subjects{cSubj},'_',' '));
    %ylim([0 600]);
    figure(retH);
    subplot(sqrt(nSubj),sqrt(nSubj),cSubj);
    hist(logtable.(subjects{cSubj}).retention_duration(logtable.(subjects{cSubj}).setsize==10));
    title(strrep(subjects{cSubj},'_',' '));
    %ylim([0 600]);
end
figure(barH);
suptitle('durations of bar, should be 250 ms');
figure(retH);
suptitle('durations of retention, should be 1250 ms');

