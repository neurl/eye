function [zscores] = zScorePupil(eye_data)
%ZSCOREPUPIL Preserve 0s as NaNs in computing z-score of pd entries.
%   Uses entire series for calculation of mean, sd. Excludes zeroes.
%   Zeroes are replaced with NaNs in the place they occurred in the series.
    pupil = eye_data.psize;
    pupil(pupil == 0) = NaN;
    mu = nanmean(pupil);
    sigma = nanstd(pupil);
    zscores = (pupil-mu)/sigma;
end

