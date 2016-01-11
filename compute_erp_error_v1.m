function erp = compute_erp_error_v1(gp, type, w, zflag, T, BL, eventName)

for sn = 1:length(gp.sub)
    data = gp.sub(sn).data;
    tm = gp.sub(sn).eye_data.tm;
    switch type
        case 'pd'
            pd = gp.sub(sn).eye_data.pd;
            for i = 1:size(pd,2)
                pd(:,i) = nanSmooth(pd(:,i), ones(w,1)/w, 'same');
            end
            pd = nanmean(pd,2);
        case 'diff_pd'
            pd = gp.sub(sn).eye_data.pd;
            for i = 1:size(pd,2)
                pd(:,i) = nanSmooth(pd(:,i), ones(w,1)/w, 'same');
            end
            pd = diff(pd')';
        case 'diff_pd_abs'
            pd = gp.sub(sn).eye_data.pd;
            for i = 1:size(pd,2)
                pd(:,i) = nanSmooth(pd(:,i), ones(w,1)/w, 'same');
            end
            pd = abs(diff(pd')');
        case 'blink'
            pd = gp.sub(sn).eye_data.blink;
            for i = 1:size(pd,2)
                pd(:,i) = nanSmooth(pd(:,i), ones(w,1)/w, 'same');
            end
    end
    if zflag
        pd = (pd - nanmean(pd))/nanstd(pd);
    end
    mess = gp.sub(sn).eye_mess;
    
    sampling_rate = ceil(1/median(diff(gp.sub(1).eye_data.tm)));
    
    % focus on first 500 (if we have 500)
    L = min(500, length(data));
    data = data(1:L);
    tID = vertcat(data.trialID);
    [tID, idx] = sort(tID(:,2));
    dt = data(idx);
    d1 = dt(1:2:end);
    d2 = dt(2:2:end);
    ag = [d1.choice] == [d2.choice];
    
    dt = [d1 d2];
    ag = [ag ag];
    for i = 1:length(dt)
        dt(i).ag = ag(i);
    end
    [~,ind] = sort([dt.trialNum]);
    data2 = dt(ind);
    
    
    
    if iscell(eventName)
        ind = strcmp(mess.str, eventName{1});
        for i = 2:length(eventName)
            ind = strcmp(mess.str, eventName{i}) | ind;
        end
    else
        ind = strcmp(mess.str, eventName);
    end
    %ind = strcmp(mess.str, 'RESP_RIGHT') | strcmp(mess.str, 'RESP_LEFT');
    et = mess.tm(ind);
    et = et(1:L);
    
    for i = 1:2
        ind = [data2.correct] == (i-1);
        iEvent = compute_iEvent(et(ind), tm);
        erp{i}(sn) = compute_erp(iEvent, pd, T, BL, 'hi', sampling_rate);
    end
    
    
end