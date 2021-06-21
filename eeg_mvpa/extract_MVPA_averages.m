function [averages, rois]= extract_MVPA_averages(stats,mask)
% extract averages from stats produced by compute_group_MVPA, based on mask
% and produces averages(roi,condition,subject)
nConditions = numel(stats);
rois = bwlabel(mask,8);
nRois = max(unique(rois));
nSubjects = size(stats(1).indivClassOverTime,1);
averages = zeros(nRois,nConditions,nSubjects);
for cCond = 1:nConditions
    for cRoi=1:nRois
        roiIndex = rois==cRoi;
        for cSubj = 1:nSubjects
            subjData = squeeze(stats(cCond).indivClassOverTime(cSubj,:,:));
            averages(cRoi,cCond,cSubj) = mean(mean(subjData(roiIndex)));
        end
    end
end