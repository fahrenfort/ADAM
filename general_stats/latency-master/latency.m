function [res,cfgNew] = latency(cfg,avgs,sign)
% [res, cfgNew] = latency(cfg,avgs,sign)
% cfg and avgs are mandatory, cfg can be empty if sign is entered as third
%   argument
% avgs is a subjects x channels x time array, contains such an array in the 
%   field 'individual' or 'data' or as specified in cfg.datafield, or
%   contains a field 'avg' which is a cell array of structures (one
%   structure for each subject)
% sign can be indicated as cfg.sign or as the third input argument
%  (argument overrules struct field); it is either 'pos'/1 or 'neg'/-1
% all timining information is in samples or in the units of cfg.times
% optional fields for cfg are (all indices as in avgs.individual):
%   extract:     which measures to extract; see the description of the
%                output [res] below
%   aggregate:   individual, GA (= grand average), or one of two jackknife
%                estimates, namely jackMiller (Miller, Patterson,& Ulrich, 
%                1998, Psycholphysiology), or jackSmulders(Smulders, 2010, 
%                Psychophysiology)
%   fig:         if true, individual ERPs and estimates are plotted
%   subs:        subject included; indices or, if cfg.subNum is given,
%                subject numbers(default = all)
%   subNum:      list of subject numbers (same order as in avgs)
%   chans:       channel included into the average; indices or, if
%                cfg.chanName is given, channel names(default = all)
%   chanName:    list of channel names (same order as in avgs)
%   peakWin:     search window for peak detection [start end]; indices or,
%                if cfg.times is given, time points(default = full range);
%                can also be a subjects x 2 matrix for individual search
%                ranges
%   meanWin:     time range across which activity is averaged for mean amplitude
%                (default = peakWin); can also be a subjects x 2 matrix
%                for individual mean windows
%   times:       all individual time points OR start/end of the submitted data;
%                if set, temporal inputs (peakWidth, peakWin, meanWin etc.)
%                are interpreted in time units instead of sampling points and
%                latency estimates (peakLat, onset, offset, width, areaLat
%                counterLat) are returned in time units
%   peakWidth:   time around peak averaged for peak amplitude; time range = 
%                peak time +- peakWidth; default = 5 sampling points or ms))
%   cWinStart:   time point from which on the counter peak is searched
%                ('peak'  or 'peakWin' (default))
%   cWinWidth:   width of the counter-peak search interval starting from
%                cfg.cWinStart (negative numbers for counter peaks
%                preceding the actual peak; indices or times);
%   cWin:        alternative to cfg.cWinStart and cfg.cWinWidth: search
%                window for the counter peak [start end]
%   percArea:    percentage of total area for percent-area latency
%                (default = 0.5)
%   percAmp:     percentage of the amplitude (peak-to-peak if counterpeak
%                is used) for on- and offset;
%   areaWin:     determines temporal boundaries for the area calculation
%                'peakWin' (default), 'ampLat', 'fullRange', or [start end]
%   areaBase:    determines one boundary in amplitude space (the other is
%                determined by the component); 'percAmp' (default if
%                calculated) or 'zero' (default otherwise) 
%   ampLatWin:   where on- and offsets are searched: 'fullRange' (default),
%                'peakWin' or [start end]
%   cBound:      Boolean that indicats whether the counter peak is one
%                border for onset/offset
%                
% res contains a (latency) measure for each subject (is a struct if 
%   several measures are requested); latencies are indices or, if cfg.time
%   is given, time points
%   mean:        mean over cfg.meanWin;
%   peakLat:     latency of the detected peak, which can also be at one border
%   onset:       time point where the ERP has reached cfg.percAmp
%   offset:      time point where the ERP has fallen back to cfg.percAmp
%   width:       difference between on- and offset
%   areaLat:     percent-area latency
%   peakAmp:     mean amplitude in the window peakLat+-peakWidth;
%   peak2peak:   difference in peak amplitude between searched peak and
%                counter peak
%   baseline:    amplitude cfg.percAmp% of the peak amplitude or in between 
%                the two peaks if counter peak is used (this was referred
%                to as 'percAmp' in ealier versions)
%   area:        sum of values above/below the cfg.areaBase
%   foundLocal:  Boolean that indicates whether a local peak was found;
%                if false, peakLat is one boundary of the search interval
%   foundOn:     Boolean that indicates whether start of the component
%                (ERP crossed the baseline) was found
%   foundOff:    Boolean that indicates whether end of the component was found
%   foundOff:    Boolean that indicates whether area latency was extracted
%                successfully
%   counterLat:  same as peakLat for the counter peak
%   counterAmp:  same as peakAmp for the counter peak
%   foundLocalC: same as foundLocal for the counter peak
% cfgNew is the same as the input cfg plus default parameters for not
% indicated parameters and parameter-sensitive explanations for each field
% of res in cfgNew.ann
%
% Typical problems:
% - instead of finding a local peak - the algorithm returns a global maximum/
%   minimum and res.foundLocal is false; increasing cfg.peakWin should help
% - the component does not return to baseline before or after the peak -
%   the respective border defined by cfg.ampLatWin is selected as on- or
%   offset and res.foundOn or res.foundOff is false, respectively
%
% Originally created by Heinrich René Liesefeld in September 2014
% Last modified by Heinrich René Liesefeld in April 2019
%
% Many thanks for bug reports and helpful suggestions to
%   Philipp Ruhnau
%   Burkhard Maess
%   Johannes Fahrenfort 
% Please send further (usability) bug reports to Heinrich.Liesefeld@psy.lmu.de

if isempty(cfg)
    cfg=struct;
end
if exist('sign','var')
    %sign of the component can also be supplied as third optional input argument
    cfg.sign=sign;
end
if isstruct(avgs)
    if isfield(cfg,'datafield')
        avgs=avgs.(cfg.datafield);
    elseif isfield(avgs,'individual')
        avgs=avgs.individual;
    elseif isfield(avgs,'data')
        avgs=avgs.data;
    else
        error('avgs must either be a subject x electrode x time array or contain a field "individual" (Fieldtrip) or "data" (EEGlab) or a user-specified field (cfg.datafield), which is such an array.')
    end
elseif iscell(avgs)
    avgsOld=avgs;
    nSubs=length(avgs);
    avgs=nan([nSubs,size(avgs{1}.avg)]);
    for subi=1:nSubs
        avgs(subi,:,:) = avgsOld{subi}.avg;
    end
end
[cfg,cfgNew]=checkConfig(cfg,avgs);

if strcmp(cfg.aggregate,'GA')
    nSubs=1;
elseif strcmp(cfg.aggregate,'jackMiller')
    nSubs=length(cfg.subs)+1;
elseif ismember(cfg.aggregate,{'individual','jackSmulders'})
    nSubs=length(cfg.subs);
end


%preallocate all variables
[comp.mean]=deal(nan(nSubs,1));
if cfg.get.peak
    [comp.peakAmp,comp.peakLat]=deal(nan(nSubs,1));
    comp.foundLocal=false(nSubs,1);
end
[comp.area,comp.baseline]=deal(nan(nSubs,1));
if cfg.get.ampLat
    [comp.onset,comp.offset,comp.width]=deal(nan(nSubs,1));
    [comp.foundOn,comp.foundOff]=deal(false(nSubs,1));
end
if cfg.get.counterPeak
    [comp.counterAmp,comp.counterLat,comp.peak2peak]=deal(nan(nSubs,1));
    comp.isLocalC=false(nSubs,1);
end
if cfg.get.areaLat
    [comp.areaLat,comp.foundArea]=deal(nan(nSubs,1));
end


if size(cfg.peakWin,1)>1 && size(cfg.peakWin,2)>1
    allPeakWin=cfg.peakWin;    
    %transpose array if it is the wrong way round
    if size(allPeakWin,1)<size(allPeakWin,2)
        allPeakWin=allPeakWin';
    end
end
if size(cfg.meanWin,1)>1 && size(cfg.meanWin,2)>1
    allMeanWins=cfg.meanWin;    
    %transpose array if it is the wrong way round
    if size(allMeanWins,1)<size(allMeanWins,2)
        allMeanWins=allMeanWins';
    end
end


for subi=1:nSubs %cycle through subjects
    if strcmp(cfg.aggregate,'GA')
        subNum=1;
    elseif strcmp(cfg.aggregate,'jackMiller') && subi==nSubs
        subNum=99;
    else
        subNum=cfg.subNum(subi);
    end
    if exist('allPeakWin','var')
        cfg.peakWin=allPeakWin(subi,:);
    end
    switch cfg.aggregate
        case 'GA'
            ERP = squeeze(mean(avgs(cfg.subs,cfg.chans,:),1));
        case {'jackMiller','jackSmulders'}
            if strcmp(cfg.aggregate,'jackMiller') && subi==nSubs
                %last value is estimate from full grand average (the value
                %actually reported)
                ERP = squeeze(mean(avgs(cfg.subs,cfg.chans,:),1));
            else
                %other values are leave-one-out (used to estimate standard
                %error with the Miller approach)
                ERP = squeeze(mean(avgs(cfg.subs(cfg.subs~=subNum),cfg.chans,:),1));
            end
        case 'individual'
            ERP = squeeze(avgs(cfg.subs(subi),cfg.chans,:));
    end
    if numel(cfg.chans) > 1 %pool over several channels
        ERP=mean(ERP,1)';
    end
    
    %mean area
    if exist('allMeanTimes','var')
        comp.mean(subi) = mean(ERP(allMeanWins(subi,1):allMeanWins(subi,2)));
    else
        comp.mean(subi) = mean(ERP(cfg.meanWin(1):cfg.meanWin(2)));
    end
    peak=peakDetection(cfg,ERP);
    comp.peakLat(subi)=cfg.times(peak.lat);
    comp.foundLocal(subi)=peak.foundLocal;
    if ~peak.foundLocal
        warning('No local peak found for subject %g',subNum)
    end
    comp.peakAmp(subi)=peak.amp;
    
    if cfg.get.counterPeak
        if isfield(cfg,'cWin')
            cPeak=peakDetection(cfg,ERP,cfg.cWin);
        else
            cPeak=peakDetection(cfg,ERP,peak.lat);
        end
        comp.counterLat(subi)=cfg.times(cPeak.lat);
        comp.counterAmp(subi)=cPeak.amp;
        comp.peak2peak(subi)=peak.amp-cPeak.amp;
        comp.isLocalC(subi)=cPeak.foundLocal;
        if ~cPeak.foundLocal
            warning('No local counter-peak found for subject %g',subNum)
        end
    end
    
    
    
    %percent amplitude
    if cfg.get.counterPeak
        %calculate threshold relative to peak-to-peak amplitude when
        %counterPeak is used
        if cfg.sign==1
            comp.baseline(subi)=peak.amp-(peak.amp-cPeak.amp)*(1-cfg.percAmp);
        elseif cfg.sign==-1
            comp.baseline(subi)=peak.amp+(cPeak.amp-peak.amp)*(1-cfg.percAmp);
        end
    else
        %calculate threshold relative to component amplitude when
        %counterPeak is not used
        comp.baseline(subi)=peak.amp-peak.amp*(1-cfg.percAmp);
    end

    if cfg.get.ampLat
        if isfield(comp,'counterLat')
            ampLat=amplitudeLatency(cfg,ERP,comp.baseline(subi),peak.lat,cPeak.lat);
        else
            ampLat=amplitudeLatency(cfg,ERP,comp.baseline(subi),peak.lat);
        end
        if cfg.warnings
            if ~ampLat.foundOn
                warning('Could not find onset for subject %g. Increase cfg.percAmp or cfg.ampLatWin to solve this problem!',subNum)
            end
            if ~ampLat.foundOff
                warning('Could not find offset for subject %g. Increase cfg.percAmp or cfg.ampLatWin to solve this problem!',subNum)
            end
        end

        comp.foundOn(subi)=ampLat.foundOn;
        comp.foundOff(subi)=ampLat.foundOff;
        comp.onset(subi)=cfg.times(ampLat.start);
        comp.offset(subi)=cfg.times(ampLat.end);
        comp.width(subi)=(ampLat.end-ampLat.start)/cfg.sampRate;

    end
    if cfg.get.area
        %get sum of all values above/below the baseline and areaLat
        switch cfg.areaBase
            case 'percAmp'
                baseline=comp.baseline(subi);
            case 'zero'
                baseline=0;
        end
        if isnumeric(cfg.areaWin)
            time=cfg.areaWin;
        else
            switch cfg.areaWin
                case 'ampLat'
                    time=[ampLat.start,ampLat.end];
                case 'peakWin'
                    time=cfg.peakWin;
                case 'fullRange'
                    time=[1,length(ERP)];
                otherwise
                    error('cfg.areaWin must be either "ampLat", "peakWin", or "fullRange", it is %s',cfg.areaWin)
            end
        end
        comp.area(subi)=totalArea(cfg,ERP,baseline,time);
        if cfg.get.areaLat
            areaLat=areaLatency(cfg,ERP,comp.area(subi),baseline,time);
        end
    end
    if cfg.get.areaLat
        if areaLat.foundLat
            comp.areaLat(subi)=cfg.times(areaLat.lat);
        else
            comp.areaLat(subi)=nan;
            warning('Did not find fractional area latency for subject %g; set latency to average of the time range',subNum);
        end
        comp.foundArea(subi)=areaLat.foundLat;
    end
    
    
    
    if cfg.fig
        if subi==1
            figure('Name','Individual ERPs and extracted estimates');
        end
        ncols=ceil(nSubs/10);
        subplot(ceil(nSubs/ncols),ncols,subi);
        if isfield(cfg,'times')
            xRange=[cfg.times(1),cfg.times(end)];
        else
            xRange=[0 length(ERP)];
        end
        yRange=[min(ERP) max(ERP)];
        plot(cfg.times,ERP,'k');
        hold on
        plot(comp.peakLat(subi),comp.peakAmp(subi),'o');
        if cfg.get.counterPeak
            plot(comp.counterLat(subi),comp.counterAmp(subi),'ro');
        end
        if cfg.get.area
            plot(xRange,[baseline baseline]);
        else
            plot(xRange,[comp.baseline(subi) comp.baseline(subi)]);
        end
        if cfg.get.ampLat
            plot([comp.onset(subi) comp.onset(subi)],yRange);
            plot([comp.offset(subi) comp.offset(subi)],yRange);
        end
        if cfg.get.areaLat
            plot([comp.areaLat(subi) comp.areaLat(subi)],yRange,'g');
        end
        plot(cfg.times([cfg.peakWin(1) cfg.peakWin(1)]),yRange,'k--');
        plot(cfg.times([cfg.peakWin(2) cfg.peakWin(2)]),yRange,'k--');
        axis('tight')
        if ismember(cfg.aggregate,{'jackMiller','jackSmulders'})
            if strcmp(cfg.aggregate,'jackMiller') && subi==nSubs
                title('Average');
            else
                title(sprintf('left out subject #%g',subNum));
            end
        else
            title(sprintf('subject #%g',subNum));
        end
    end
    
end %end of loop over subjects

if exist('allPeakWin','var')
    cfg.peakWin=allPeakWin;
end

if mean(comp.foundLocal) < 0.5
    %check the proportion of local peaks found
    fprintf('Attention: In less than half of the cases a local peak was found.\n');
    fprintf('Consider increasing the search interval.\n');
end
if isfield(comp,'isLocalC') && mean(comp.isLocalC) < 0.5
    %check the proportion of local counter peaks found
    fprintf('Attention: In less than half of the cases a local counter peak was found.\n');
    fprintf('Consider increasing the search interval.\n');
end
if isfield(comp,'foundOn') && mean(comp.foundOn) <0.5
    fprintf('Attention: In less than half of the cases the onset of the component was found.\n');
    fprintf('Consider increasing the on-/offset search interval (cfg.ampLatWin), changing cfg.percAmp, and/or using the counter-peak method.\n');
end
if isfield(comp,'foundOff') && mean(comp.foundOff) <0.5
    fprintf('Attention: In less than half of the cases the offset of the component was found.\n');
    fprintf('Consider increasing the on-/offset search interval (cfg.ampLatWin), changing cfg.percAmp, and/or using the counter-peak method.\n');
end
if isfield(cfg,'chanName')
    cfg.chans=cfg.chanName(cfg.chans);
elseif any(ischar(cfg.chans)) || any(mod(cfg.chans,1))
    error('If cfg.chanName is not set, cfg.chans must contain indices, not channel names.')
end
if iscell(cfg.extract)
    for outi=1:length(cfg.extract)
        if strcmp(cfg.aggregate,'jackSmulders')
            %oi = n*mean(J) - (n - 1)*ji, see Smulders (2010)
            for subi=1:nSubs
                res.(cfg.extract{outi})(subi)=nSubs*mean(comp.(cfg.extract{outi}))-(nSubs-1)*comp.(cfg.extract{outi})(subi);
            end
        else
            res.(cfg.extract{outi})=comp.(cfg.extract{outi});
        end
    end
elseif strcmp(cfg.aggregate,'jackSmulders')
    %oi = n*mean(J) - (n - 1)*ji, see Smulders (2010)
    for subi=1:nSubs
        res(subi)=nSubs*mean(comp.(cfg.extract))-(nSubs-1)*comp.(cfg.extract)(subi);
    end
    
else
    res=comp.(cfg.extract);
end
end

function [cfg,cfgNew]=checkConfig(cfg,avgs)
%cfg.sign must be indicated or function returns an error
%avgs is necessary when cfg.subs or cfg.chans or cfg.peakWin is
%missing (to determine the maximum extend in each of these dimensions)

cfgNew=[];
if ~isfield(cfg,'sign')
    error('Indicating the sign of the component (cfg.sign = ''pos'' or 1 or ''neg'' or -1) is mandatory.')
end
if ischar(cfg.sign)
    switch cfg.sign
        case 'pos'
            cfg.sign=1;
        case 'neg'
            cfg.sign=-1;
    end
end

if ~isfield(cfg,'warnings'),cfg.warnings=1;end
if ~isfield(cfg,'subs') || strcmp(cfg.subs,'all')
    cfg.subs=1:size(avgs,1);
end
if isfield(cfg,'subs')
    if all(ismember(cfg.subs,[1,0]))
        cfg.subs = logical(cfg.subs);
    end
    if islogical(cfg.subs)
        %cfg.subs can also be set as a filter
        cfg.subs=find(cfg.subs);
    end
elseif isfield(cfg,'subNum')
    allSubs=1:length(cfg.subNum);
    cfg.subs=allSubs(ismember(cfg.subNum,cfg.subs));
end
if ~isfield(cfg,'subNum')
    cfg.subNum=cfg.subs;
end

if ~isfield(cfg,'chans')
    if cfg.warnings
        warning('No channels (cfg.chans) indicated; will average across all of them.')
    end
    cfg.chans=1:size(avgs,2);
elseif isfield(cfg,'chanName') && (iscell(cfg.chans) || ischar(cfg.chans))
    allChans=1:length(cfg.chanName);
    cfg.chans=allChans(ismember(cfg.chanName,cfg.chans));    
end
if ~isfield(cfg,'peakWin')
    if cfg.warnings
        warning('No time range (cfg.peakWin) indicated; will use full time range.');
    end
    cfg.peakWin=[1,size(avgs,3)];
end

if isfield(cfg,'times')
    %determine format of time
    if ~isfield(cfg,'timeFormat')
        oneSample=cfg.times(2)-cfg.times(1);
        if oneSample > 0.1 %likely in ms (if this was in s, sampling rate would be < 10 samples/s)
            cfg.timeFormat='ms';
            fprintf('All times are interpreted as milliseconds; set cfg.timeFormat=''s'' if this is incorrect.\n')
        else %likely in s (if this was in ms, sampling rate would be >= 10 samples/ms)
            cfg.timeFormat='s';
            fprintf('All times are interpreted as seconds; set cfg.timeFormat=''ms'' if this is incorrect.\n')            
        end
    end
    
    %check cfg.times for consistency
    if length(cfg.times)==1
        error('cfg.times must include at least two values (start and end times).')
    elseif length(cfg.times)==3
        cfg.times=cfg.times(1:2);
        if cfg.warnings
            warning('Sampling rate is now calculated from the data; third value in cfg.times is ignored.');
        end
    end
    if length(cfg.times)==2
        cfg.times=linspace(cfg.times(1),cfg.times(2),size(avgs,3));
    elseif length(cfg.times) ~= size(avgs,3)
        error('Length of cfg.times must match the number of samples in avgs!')
    end
    if any(diff(cfg.times)-mean(diff(cfg.times))~=0)
        %non-constant sampling rate
        switch cfg.timeFormat
            %deviations in sub-microsecond range do not matter
            case 'ms'
                tolerance=1e-3;
            case 's'
                tolerance=1e-6;
        end
        if all(abs(diff(cfg.times)-mean(diff(cfg.times)))<tolerance)
            cfg.times=round(cfg.times/tolerance)*tolerance;
            if cfg.warnings
                warning('Sligthly rounded cfg.times, to achieve constant sampling rate.')
            end
        else
            error('Apparently non-constant sampling rate - check cfg.times!')
        end
    end
    
    %determine sampling rate
    if isfield(cfg,'sampRate') && cfg.warnings
        warning('Sampling rate is now calculated from the data; cfg.sampRate is ignored.');
    end
    cfg.sampRate=1/mean(diff(cfg.times));
    fprintf('Calculated sampling rate: %g samples/%s\n',cfg.sampRate,cfg.timeFormat);
else
    fprintf('All times are interpreted as sampling points; fill cfg.times to change this.\n')
    cfg.timeFormat = 'samples';
    cfg.times=1:size(avgs,3);
    cfg.sampRate=1;
    tolerance=0;
end

if ~isfield(cfg,'meanWin')
    if isfield(cfg,'meanTime')
        warning('cfg.meanTime is now called cfg.meanWin for consistency.')
        cfg.meanWin=cfg.meanTime;
    else
        cfg.meanWin=cfg.peakWin;
    end
end

%check and/or determine cfg.peakWidth
if isfield(cfg,'peakWidth')
    cfg.peakWidth=round(cfg.peakWidth*cfg.sampRate);
    cfgNew.peakWidth=cfg.peakWidth/cfg.sampRate;
    if cfg.warnings && abs(cfgNew.peakWidth-cfg.peakWidth)>tolerance
        warning('Rounded peakWidth (%g %s) to the next possible value: %g %s',...
            cfg.peakWidth,cfg.timeFormat,cfgNew.peakWidth,cfg.timeFormat);
    end
else
    %set default peak width: ~5ms or 5 samples
    switch cfg.timeFormat
        case 's'
            peakWidthS=0.005; %default to ~5 ms
            cfg.peakWidth=round(peakWidthS*cfg.sampRate);
            peakWidthS=cfg.peakWidth/cfg.sampRate;
            if cfg.warnings
                warning('No peak width (cfg.peakWidth) indicated; peakWidth set to %g s.',peakWidthS)
            end
        case 'ms'
            peakWidthMs=5; %default to ~5 ms
            cfg.peakWidth=round(peakWidthMs*cfg.sampRate);
            peakWidthMs=cfg.peakWidth/cfg.sampRate;
            if cfg.warnings
                warning('No peak width (cfg.peakWidth) indicated; peakWidth set to %g ms.',peakWidthMs)
            end
        case 'samples' %default to 5 sampling points
            cfg.peakWidth=5;
            if cfg.warnings
                warning('No peak width (cfg.peakWidth) indicated; peakWidth set to %g asmpling points.',cfg.peakWidth)
            end
    end
end

%handle all the timing information
checkFields={'cWinWidth','peakWin','meanWin','cWin','areaWin','ampLatWin'};
if ismember(cfg.timeFormat,{'s','ms'})
    %transform to sampling points
    for i=1:length(checkFields)
        if isfield(cfg,checkFields{i}) && isnumeric(cfg.(checkFields{i}))
            if isempty(regexp(checkFields{i},'Width','once'))
                [samplesOut,delta] = deal(zeros(size(cfg.(checkFields{i}))));
                maxDelta      = max(diff(cfg.times)); %because it is not always perfectly equidistant
                for j=1:length(cfg.(checkFields{i}))
                    [delta(j),samplesOut(j)] = min(abs(cfg.times-cfg.(checkFields{i})(j)));
                    if delta(j) > maxDelta
                        error('Time point %g %s (%s) not found in time range [%g - %g %s]; maybe something is wrong about (the units of) cfg.times',...
                            cfg.(checkFields{i})(j),cfg.timeFormat,checkFields{i},cfg.times(1),cfg.times(end),cfg.timeFormat);
                    end
                end
                if cfg.warnings && any(abs(cfg.(checkFields{i})-cfg.times(samplesOut))>tolerance)
                    if length(cfg.(checkFields{i}))==2
                        warning('Rounded %s (%g-%g %s) to the next possible values: %g - %g %s',...
                            checkFields{i},cfg.(checkFields{i}),cfg.timeFormat,cfgNew.(checkFields{i}),cfg.timeFormat);
                    else
                        warning('Rounded %s (%g %s) to the next possible value: %g %s',...
                            checkFields{i},cfg.(checkFields{i}),cfg.timeFormat,cfgNew.(checkFields{i}),cfg.timeFormat);
                    end
                end            
                cfgNew.(checkFields{i})=cfg.times(samplesOut);
            else
                samplesOut=round(cfg.(checkFields{i})*cfg.sampRate);
                cfgNew.(checkFields{i})=samplesOut/cfg.sampRate;
                if cfg.warnings && abs(cfgNew.(checkFields{i})-cfg.(checkFields{i}))>tolerance
                    warning('Rounded %s (%g %s) to the next possible value: %g %s',...
                        checkFields{i},cfg.(checkFields{i}),cfg.timeFormat,cfgNew.(checkFields{i}),cfg.timeFormat);
                end
            end
            cfg.(checkFields{i})=samplesOut;
        end
    end
else
    %check whether all timing information is indicated in sampling points
    msg='If cfg.times is not set, all timing information is interpreted as sampling points and must consequently be in integers. Please check: ';
    throwError=false;
    for i=1:length(checkFields)
        if isfield(cfg,checkFields{i}) && isnumeric(cfg.(checkFields{i}))
            if any(mod(cfg.(checkFields{i}),1)) || (isempty(regexp(checkFields{i},'Width','once')) && (any(cfg.(checkFields{i})<1) || any(cfg.(checkFields{i})>size(avgs,3))))
                msg=[msg,' ',checkFields{i}];
                throwError=true;
            end
        end
    end
    if throwError
        error(msg)
    end
end

%check cfg.extract (requested output measures)
extractable={'peakLat','onset','offset','areaLat','mean','peakAmp','area','width','peak2peak','baseline'};
if ~isfield(cfg,'extract')
    cfg.extract=extractable;
    if cfg.warnings
        warning('No output measures chosen; will extract all possible output measures:')
        fprintf('%s, ',cfg.extract{:});
        fprintf('\n');
    end
end
if strcmp(cfg.extract,'all')
    cfg.extract=extractable;
end
if ismember('percAmp',cfg.extract)
    if iscell(cfg.extract)
        cfg.extract{ismember(cfg.extract,'percAmp')}='baseline';
    else
        cfg.extract='baseline';
    end
    warning('''percAmp'' is now called ''baseline''; call ''help latency'' or see cfgNew.ann.baseline for details')
end
if any(~ismember(cfg.extract,extractable))
    fprintf('\nAvailable output measures:\n');
    fprintf('%s\n',extractable{:});
    fprintf('call ''help latency'' for details\n');
    if iscell(cfg.extract)
        error('Cannot extract ''%s''\n',cfg.extract{~ismember(cfg.extract,extractable)});
    else
        error('Cannot extract ''%s''',cfg.extract)
    end
end

%determine which calculations to perform
cfg.get.peak=true; %needs nothing in addition

if ismember('areaLat',cfg.extract)
    if ~isfield(cfg,'areaWin')
        if isfield(cfg,'percAmp') && (~isfield(cfg,'areaBase') || ~strcmp(cfg.areaBase,'zero'))
            cfg.areaWin='ampLat';
        else
            cfg.areaWin='peakWin';
        end
    end
    if ~isfield(cfg,'percArea'),cfg.percArea=.5;end
    cfg.get.areaLat=true;
else
    cfg.get.areaLat=false;
end

if any(ismember({'onset','offset','width'},cfg.extract)) || (cfg.get.areaLat && strcmp(cfg.areaWin,'ampLat'))
    cfg.get.ampLat=true; %needs only peak
    if ~isfield(cfg,'ampLatWin'),cfg.ampLatWin='fullRange';end
else
    cfg.get.ampLat=false;
end

if any(ismember({'areaLat','area'},cfg.extract))
    cfg.get.area=true;
else
    cfg.get.area=false;
end
%searches for the point where the set percentage  of the total area is
%reached and therefore critically depends on area
if ~isfield(cfg,'fig'),cfg.fig=false;end
if ~isfield(cfg,'aggregate')
    if isfield(cfg,'aggregation')
        cfg.aggregate=cfg.aggregation;
        cfg=rmfield(cfg,'aggregation');
    else
        cfg.aggregate='individual';
    end
end
if isfield(cfg,'cWinWidth')
    if ~isfield(cfg,'cWinStart')
        cfg.cWinStart='peakWin';
    end
    %if set, percAmp determines the peak-to-peak threshold instead of the
    %x-axis-to-peak threshold
    if ~isfield(cfg,'cBound'),cfg.cBound=true;end
    %search window for releative latency is bounded on one side by the counter peak
    if cfg.cBound,cfg.get.counterPeak=1;end
elseif isfield(cfg,'cWin')
    cfg.get.counterPeak=1;
else
    cfg.get.counterPeak=0;
end

if isfield(cfg,'percAmp')
    cfg.get.newBaseline=1;
    %the area is bound in amplitude by percAmp of the peak amplitude or of
    %the peak-to-peak amplitude (when counterPeak is calculated)
    if ~isfield(cfg,'areaBase'),cfg.areaBase='percAmp';end
else
    %otherwise area uses the typical baseline
    cfg.get.newBaseline=0;
    if any(ismember({'onset','offset','width'},cfg.extract))
        cfg.percAmp=0.5; %because 0 does not make sense (will be caught by noise)
    else
        cfg.percAmp=0; %makes sense for areaLat
    end
    if ~isfield(cfg,'areaBase'),cfg.areaBase='zero';end
end
if isnumeric(cfg.areaBase) && cfg.areaBase==0
    cfg.areaBase='zero';
end


if cfg.sign==1
    direction='posi';
elseif cfg.sign==-1
    direction='nega';
end


%annotations that explain in some detail what each measure reflects
cfg.ann.mean=sprintf('Mean activity in the interval %g-%g %s',...
    cfg.times(cfg.meanWin),cfg.timeFormat);
cfg.ann.peak=sprintf('Most %stive peak in the search interval %g-%g %s',...
    direction,cfg.times(cfg.peakWin),cfg.timeFormat);
if cfg.get.counterPeak
    if isfield(cfg,'cWin')
        cfg.ann.counterPeak=sprintf('Strongest peak of opposite polarity in the time window %g - %g %s',...
            cfg.times(cfg.cWin),cfg.timeFormat);
    else
        if cfg.cWinWidth<0 && strcmp(cfg.cWinStart,'peak')
            cfg.ann.counterPeak=sprintf('Strongest peak of opposite polarity in between %g samples before the actual peak and the actual peak',-cfg.cWinWidth);
        elseif cfg.cWinWidth>0 && strcmp(cfg.cWinStart,'peak')
            cfg.ann.counterPeak=sprintf('Strongest peak of opposite polartity in between the actual peak and %g samples after the actual peak',cfg.cWinWidth);
        elseif cfg.cWinWidth<0 && strcmp(cfg.cWinStart,'peakWin')
            cfg.ann.counterPeak=sprintf('Strongest peak of opposite polarity in between %g samples before the start of the search window and the actual peak',-cfg.cWinWidth);
        elseif cfg.cWinWidth>0 && strcmp(cfg.cWinStart,'peakWin')
            cfg.ann.counterPeak=sprintf('Strongest peak of opposite polartity in between the actual peak and %g samples after the end of the search window',cfg.cWinWidth);
        end
    end
end
if cfg.get.newBaseline
    if cfg.get.counterPeak
        cfg.ann.baseline=sprintf('Activity %g%% of the peak-to-counter-peak distance away from the peak\n(lower values are closer to the peak)',...
            cfg.percAmp*100);
    else
        cfg.ann.baseline=sprintf('Activity %g%% of the peak amplitude away from the peak amplitude\n(lower values are closer to the peak)',...
            cfg.percAmp*100);
    end
    cfg.ann.onset=sprintf('Point before the peak, where activity crosses the baseline\n(see cfg.ann.baseline)');
    cfg.ann.offset=sprintf('Point after the peak, where activity crosses the baseline\n(see cfg.ann.baseline)');
    
    cfg.ann.area=sprintf('Cumulative sum over all %stive samples in between on- and offset (see cfg.ann.onset),\nthe ERP and the %g%% threshold (see cfg.ann.baseline)',...
        direction,cfg.percAmp*100);
else
    cfg.ann.area=sprintf('Cumulative sum over all %stive samples in between on- and offset (see cfg.ann.onset),\nthe ERP and the pre-stimulus baseline activity',...
        direction);
end
if cfg.get.areaLat
    cfg.ann.areaLat=sprintf('Point in time where %g%% of the total area (see cfg.ann.area) is reached',cfg.percArea*100);
end

fields=fieldnames(cfg);
for fieldi=1:length(fields)
    if ~isfield(cfgNew,fields{fieldi})
        cfgNew.(fields{fieldi})=cfg.(fields{fieldi});
    end
end

end

function peak=peakDetection(cfg,ERP,peakLat)
%will search for counterpeak if peakLat is set

if nargin<3
    peakWin=cfg.peakWin;
    direction=cfg.sign;
else
    %counter peak has the opposite polarity
    direction = cfg.sign*-1;
    if length(peakLat)==1
        switch cfg.cWinStart
            case 'peak'
                if cfg.cWinWidth < 0 %preceding counter peak
                    peakWin=[peakLat+cfg.cWinWidth,peakLat];
                elseif cfg.cWinWidth >= 0 %following counter peak
                    peakWin=[peakLat,peakLat+cfg.cWinWidth];
                end
            case 'peakWin'
                if cfg.cWinWidth < 0 %preceding counter peak
                    peakWin=[cfg.peakWin(1)+cfg.cWinWidth,peakLat];
                elseif cfg.cWinWidth >= 0 %following counter peak
                    peakWin=[peakLat,cfg.peakWin(2)+cfg.cWinWidth];
                end
        end
    elseif length(peakLat)==2
        peakWin=peakLat;
    end
end

%make sure not to get a peak that lies directly at one of the ends of the
%ERP so that amplitude can always be calculated
peakWin(1)=max(peakWin(1),cfg.peakWidth+1);
peakWin(2)=min(peakWin(2),length(ERP)-cfg.peakWidth);
    
ERPcut = ERP(peakWin(1):peakWin(2));

if direction==-1
    ERPcut=ERPcut*-1;
    %if a negative peak is searched, search for local positive peak on
    %the reversed  polarity ERP
end
if length(ERPcut)>2
    [dummy,lats] = findpeaks(ERPcut,'SORTSTR', 'descend');
    if numel(lats)<1
        %if no local peak is found (this is what findpeaks does), simply take the maximum
        [dummy,lat]=max(ERPcut);
        %track whether a local peak was found or not
        peak.foundLocal= false;
    else %take the strongest local peak
        lat = lats(1);
        peak.foundLocal = true;
    end
else
    if cfg.warnings
        warning('Search range for peak is too narrow (less than 3 samples)!')
    end
    [dummy,lat]=max(ERPcut);
    peak.foundLocal= false;
end

%convert from relative to search range TO relative to episode onset
peak.lat = peakWin(1)+lat-1;

peak.amp = mean(ERP(peak.lat-cfg.peakWidth:peak.lat+cfg.peakWidth));
%use original data instead of cutERP, because the peakWin might overlap
%with start/end of cutERP
end

function ampLat=amplitudeLatency(cfg,ERP,baseline,peakLat,counterLat)
%ampLat.start: where ERP reaches cfg.percAmp before the peak
%ampLat.end:   where ERP reaches cfg.percAmp after the peak


%determine the minimal and maximal values for ampLat
%in case of peakWin, one side is bounded by the counter peak (if there is one)
if ischar(cfg.ampLatWin)
    switch cfg.ampLatWin
        case 'peakWin'
            time=cfg.peakWin;
        case 'fullRange'
            time=[1,length(ERP)];
        otherwise
            error('cfg.ampLatWin must either be "peakWin", "fullRange" or start and end time indices; it is %s',cfg.ampLatWin)            
    end
    if exist('counterLat','var') && cfg.cBound
        %counter peak is one border of the interval
        if counterLat<peakLat
            time(1)=counterLat;
        elseif counterLat>=peakLat
            time(2)=counterLat;
        end
        
    end
elseif isnumeric(cfg.ampLatWin) && numel(cfg.ampLatWin)==2;
    time=cfg.ampLatWin;
else
    error('cfg.ampLatWin must either be "peakWin", "fullRange" or start and end time indices')
end

%make sure to avoid samples that lie directly at one of the ends of the
%ERP so that we certainly can calculate an amplitude
time(1)=max(time(1),cfg.peakWidth+1);
time(2)=min(time(2),length(ERP)-cfg.peakWidth);


%first search for early crossing of baseline by going backwards from
%peak
ampLat=struct;
for sample=peakLat:-1:time(1)+cfg.peakWidth
    %go through the search window and calculate mean area (peakWidth
    %broad) and each time check whether baseline has been reached
    currentAmp = mean(ERP(sample-cfg.peakWidth:sample+cfg.peakWidth));
    if (cfg.sign==1 && currentAmp<=baseline) || (cfg.sign==-1 && currentAmp>=baseline)
        ampLat.start=sample;
        break;
    end
    
end

%now search for the later crossing of the threshold
for sample=peakLat:time(2)-cfg.peakWidth
    currentAmp = mean(ERP(sample-cfg.peakWidth:sample+cfg.peakWidth));
    if (cfg.sign==1 && currentAmp<=baseline) || (cfg.sign==-1 && currentAmp>=baseline)
        ampLat.end=sample;
        break;
    end
end

%When ERP did not cross the baseline in the search interval,
%percent-amplitude latency is set to the borders of the interval
if isfield(ampLat,'start')
    ampLat.foundOn=true;
else
    ampLat.foundOn=false;
    ampLat.start=time(1);
end
if isfield(ampLat,'end')
    ampLat.foundOff=true;
else
    ampLat.foundOff=false;
    ampLat.end=time(2);
end

end %end of sub-function ampLat

function area=totalArea(cfg,ERP,baseline,time)
cutERP = (ERP(time(1):time(2))-baseline)*cfg.sign;
area = sum(cutERP(cutERP>0))*cfg.sign; % area is negative if the component is negative
end %end of sub-function totalArea

function areaLat=areaLatency(cfg,ERP,area,baseline,time)
%time in samples (taken care of by checkConfig)
cutERP=(ERP(time(1):time(2))-baseline)*cfg.sign; % ERP & area need to be flipped if they are negative
sum_tempERP = 0;
for sample=1:size(cutERP)
    %as long as we have not yet reached the quantile area
    %specified by percArea, go on to add the value at the
    %current sample to the sum; else break (because
    %percArea-%-area latency is reached
    if cutERP(sample) > 0
        if sum_tempERP < area*cfg.sign*cfg.percArea
            sum_tempERP = sum_tempERP+cutERP(sample);
        else
            areaLat.lat = sample+time(1)-1;
            break;
        end
    end
end
% does not work: areaLat.lat = find(cumsum(cutERP(cutERP>0)) >= area * cfg.sign * cfg.percArea,1) + time(1)-1;
if exist('areaLat','var')
   areaLat.foundLat=true;
else
    areaLat.lat = NaN;
    areaLat.foundLat=false;
end
end %end of sub-function areaLat



