function iEvent = compute_iEvent(t1, t2)
% find indices of events in pupil data
%?t1 = message times
%?t2 = pupil time series
count = 1;
iEvent = nan(1,length(t1));
dum_old = 0;
for i = 1:length(t1)
    %return the greatest value of t2 less than t1(i)
    dum = find(t1(i) > t2(dum_old+1:end), 1, 'last') + dum_old;
    if ~isempty(dum)
        iEvent(count) = dum;
        count = count + 1;
    end
    dum_old = dum;
end
iEvent(isnan(iEvent)) = [];