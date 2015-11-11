function [pdp, idx] = remove_pd_outliers(pd, thresh)


% remove top thresh and bottom thresh fraction of pds
[s, ind] = sort(pd);
%exclude values that are already nan
ind(isnan(s))=[];
s(isnan(s))=[];
t = [1:length(s)]/sum(~isnan(s));
sp = s;
sp((t<thresh) | (t>(1-thresh))) = nan;
pdp(ind) = sp;
idx = isnan(pdp);