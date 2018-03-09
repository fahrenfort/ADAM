function [ subjnames_or_datastruct, inclindex, exclindex ] = select_subjects(subjnames_or_datastruct,inclsubj,invert)
% function [ subjnames_or_datastruct, inclindex, exclindex ] = select_subjects(subjnames_or_datastruct,inclsubj,invert)
% function to select subjects specified by inclsubj
% subjnames_or_datastruct is either a weightstruct (function can later
% easily be expanded to include other datastructs such as statstructs) or a
% cell array of subject names
% inclsubj either contains the index numbers of the subjects to include, a
% cell array of the subject names to be included, or a string pattern to
% search for
% when invert is true, the selection is inverted (removing subjects in
% inclsubj rather than including them)
% J.J.Fahrenfort, VU, 2016

if nargin<3
    invert = false;
end
if nargin<2
    error('requires two inputs: subjnames_or_datastruct and inclsubj');
end

if iscell(subjnames_or_datastruct) % selecting names from cell array containg subjects
    nSubj = numel(subjnames_or_datastruct);
    if ~isempty(inclsubj)
        if ischar(inclsubj)
            inclsubj = {inclsubj};
        end
        if iscell(inclsubj)
            inclsubjIndex = [];
            for cPattern = 1:numel(inclsubj)
                subjIndex = strfind(subjnames_or_datastruct,inclsubj{cPattern});
                inclsubjIndex = [inclsubjIndex find(not(cellfun('isempty', subjIndex)))];
            end
            inclsubj = unique(inclsubjIndex);
        else
            inclsubj = find(ismember(1:nSubj,inclsubj));       
        end
        exclsubj = find(~ismember(1:nSubj,inclsubj));
        if invert
            [inclsubj, exclsubj] = swapvars(inclsubj, exclsubj);
        end
        disp(['including subject(s): ' strjoin(subjnames_or_datastruct(inclsubj),', ')]);
        disp(['excluding subject(s): ' strjoin(subjnames_or_datastruct(exclsubj),', ')]);
        subjnames_or_datastruct = subjnames_or_datastruct(inclsubj);
        inclindex = inclsubj;
        exclindex = exclsubj;
    else
        disp('no subjects are removed/selected');
        inclindex = 1:nSubj;
        exclindex = [];
    end
else % or removing data from datastruct
    for cStruct = 1:numel(subjnames_or_datastruct)
        substr = subjnames_or_datastruct(cStruct);
        fld = fieldnames(substr)';
        if isfield(substr,'filenames')
            nSubj = numel(substr.filenames);
            [ subjnames, inclindex, exclindex ] = select_subjects(substr.filenames,inclsubj,invert);
        else
            nSubj = size(substr.(fld{find(strncmp(fld,'indiv',5),1)}),1);
            if iscell(inclsubj) || ischar(inclsubj)
                error('no filenames field present, cannot select subjects based on name or pattern');
            end
            inclindex = inclsubj;
            exclindex = find(~ismember(1:nSubj,inclindex));
        end
        if numel(exclindex) == 0 || numel(exclindex) == nSubj
            disp('WARNING: selection unsuccessful (no or all subjects were selected), keeping all subjects by default.');
            inclindex = 1:nSubj;
            exclindex = [];
        else
            if isfield(substr,'filenames')
                disp(['including subjects: ' strjoin(subjnames,', ')]);
                substr.filenames = subjnames;
            else
                disp(['including subjects: ' regexprep(num2str(exclindex),' +',', ')]);
            end
        end
        for field = fld
            fld = field{1};
            if strcmp(fld,'indivCTF')
                disp(['limiting subjects in ' fld ' and recomputing CTF and semCTF']);
                dat = substr.(fld);
                dat = dat(inclindex,:,:,:,:);
                substr.(fld) = dat;
                substr.CTF = squeeze(mean(dat,1));
                substr.semCTF = squeeze(std(dat)/sqrt(numel(inclindex)));
            elseif strcmp(fld,'indivCTFpercond')
                disp(['limiting subjects in  ' fld ' and recomputing CTFpercond and semCTFpercond']);
                for cCond = 1:numel(substr.(fld))
                    dat = substr.(fld){cCond};
                    dat = dat(inclindex,:,:,:,:);
                    substr.(fld){cCond} = dat;
                    substr.CTFpercond{cCond} = squeeze(mean(dat,1));
                    substr.semCTFpercond{cCond} = squeeze(std(dat)/sqrt(numel(inclindex)));
                end
            elseif strcmp(fld,'indivWeights')
                disp(['limiting subjects in  ' fld ' and recomputing avWeights']);
                dat = substr.(fld);
                dat = dat(inclindex,:,:,:,:);
                substr.(fld) = dat;
                substr.avWeights = squeeze(mean(dat,1));
            elseif strcmp(fld,'indivClassOverTime')
                disp(['limiting subjects in  ' fld ' and recomputing ClassOverTime and StdError']);
                dat = substr.(fld);
                dat = dat(inclindex,:,:,:,:);
                substr.(fld) = dat;
                substr.ClassOverTime = squeeze(mean(dat,1));
                substr.StdError = squeeze(std(dat)/sqrt(numel(inclindex)));
                substr.pVals = [];
            elseif strcmp(fld,'indivClassAv')
                dat = substr.(fld);
                dat = dat(inclindex);
                substr.(fld) = dat;
            end
        end
        subjnames_or_datastruct(cStruct) = substr;
    end
    disp('done!');
end