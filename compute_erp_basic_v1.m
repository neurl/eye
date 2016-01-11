function erp = compute_erp_basic_v1(gp, type, w, zflag, T, BL, eventName)

for sn = 1:length(gp.sub)
    sampling_rate = ceil(1/median(diff(gp.sub(1).eye_data.tm)));
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
    
    % clicks on
    
    %T = [-2 5];
    %BL = [-1 0];
    
    ind = strcmp(mess.str, eventName);
    %ind = strcmp(mess.str, 'RESP_LEFT') | strcmp(mess.str, 'RESP_RIGHT');
    %sm(sn) = sum(ind);
    %ind = strcmp(mess.str, 'OUT_CORRECT') | strcmp(mess.str, 'OUT_WRONG');
    et = mess.tm(ind);
    iEvent = compute_iEvent(et, tm);
    erp(sn) = compute_erp(iEvent, pd, T, BL, 'hi', sampling_rate);
end