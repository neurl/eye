function pd = clean_pupil(eye_data_struct)

% THIS MIGHT NEED SOME WORK TO MAKE IT REALLY ROBUST
% reject trials when fixation is broken
ed = eye_data_struct;
pd = ed.psize;

% turn dropouts into nans (mean diameter)
pd(pd==0) = nan;
% turn dropouts into nans (individual eye)
ed.Lpsize(ed.Lpsize==0) = nan;
ed.Rpsize(ed.Rpsize==0) = nan;

% if either eye is missing remove other eye too
ind = (isnan(ed.Lpsize)|isnan(ed.Rpsize));
pd(ind) = nan;

% remove points where two pupil diameters are out of whack
thresh = 0.02;
[~, idx] = remove_pd_outliers(ed.Lpsize-ed.Rpsize, thresh);
pd(idx) = nan;
% pdp = remove_pd_outliers(pd, thresh)

end