function erp = compute_erp(iEvent, pd, T, BL, name, sampling_rate)


if isempty(iEvent)
    D = nan(1, ...
        max((post_time-pre_time)*sampling_rate)+1);
    Db = nan(1, ...
        max((post_time-pre_time)*sampling_rate)+1);
    tEvent = nan;
    tEnd = nan;
    return;
end

if ~exist('sampling_rate')
    sampling_rate = 60;
end

pre_time = T(1);
post_time = T(2);

bl_start = BL(1);
bl_end = BL(2);

% % find indices of events in pupil data
% count = 1;
% iEvent = nan(1,length(t1));
% dum_old = 0;
% for i = 1:length(t1)
%     dum = find(t1(i) > t2(dum_old+1:end), 1, 'last') + dum_old;
%     if ~isempty(dum)
%         iEvent(count) = dum;
%         count = count + 1;
%     end
%     dum_old = dum;
% end
% iEvent(isnan(iEvent)) = [];




tEvent = sampling_rate*pre_time;
tEnd = sampling_rate*(pre_time+post_time);

iStart = iEvent + sampling_rate*pre_time;
iEnd = iEvent + sampling_rate*post_time;

iEnd(iEnd > length(pd)) = length(pd);
iStart(iStart < 1) = 1;

iBL_start = iEvent + bl_start * sampling_rate;
iBL_end = iEvent + bl_end * sampling_rate;

iBL_start(iBL_start < 1) = 1;
iBL_end(iBL_end > length(pd)) = length(pd);

% get data
D = nan(max(iEnd-iStart)+1, size(pd,2), length(iStart));
for i = 1:length(iStart)
    ind = iStart(i):iEnd(i);
    D(1:length(ind),:,i) = pd(ind,:);
    
    BL_ind = iBL_start(i):iBL_end(i);
    BL(i,:) = nanmean(pd(BL_ind,:));
    for j = 1:size(pd,2)
        Db(1:length(ind),j,i) = D(1:length(ind),j,i) - BL(i,j);
    end
end

tm = [T(1):1/sampling_rate:T(2)];
erp.D           = D;
erp.Db          = Db;
erp.y           = squeeze(nanmean(Db,3));
erp.y_nb        = squeeze(nanmean(D,3));
erp.type        = name;
erp.tEvent      = tEvent;
erp.tEnd        = tEnd;
erp.tm          = tm;
