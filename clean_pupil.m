function pd_clean = clean_pupil(pd_raw)

% THIS MIGHT NEED SOME WORK TO MAKE IT REALLY ROBUST
% reject trials when fixation is broken

pd = pd_raw;

% turn dropouts into nans
pd(pd==0) = nan;

% if either eye is missing remove other eye too
ind = (isnan(pd(:,1))|isnan(pd(:,2)));
pd(ind,:) = nan;

% remove points where two pupil diameters are out of whack
thresh = 0.02;
[~, idx] = remove_pd_outliers(pd(:,1)-pd(:,2), thresh);
pd(idx, :) = nan;
% pdp = remove_pd_outliers(pd, thresh)

% average R and L pupil diameters
pd_clean = mean(pd(:,1:2),2);
