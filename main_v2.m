%% setup
clear

maindir = ['~/Work/Projects/clicks/task_clicks/'];
datadir = [maindir 'datadir/2015_frozenNoise_v3/'];
fundir  = [maindir 'fundir/v1/'];

addpath(fundir);
cd(fundir);

defaultPlotParameters

gp = group(datadir);
gp.load;

%% clean pupil
gp.clean_pupil;

%% exclude bad subjects
clear perf dropoutRate
for sn = 1:length(gp.sub)
    perf(sn) = mean([gp.sub(sn).data.correct]);
    dropoutRate(sn) = mean(mean(isnan(gp.sub(sn).eye_data.pd)));
end
figure(1); clf;
hist(perf)

% exclude below threshold
thresh = 0.55;
ind_bad = (perf < thresh) | (dropoutRate > 0.5);
gp.sub = gp.sub(~ind_bad);

% based on eye data






%% detect blinks
for sn = 1:length(gp.sub)
    t_min = 50/1000; t_max = 500/1000;
    
    samplingRate = 1/(gp.sub(sn).eye_data.tm(2)-gp.sub(sn).eye_data.tm(1));
    pd = mean(gp.sub(sn).eye_data.pd_raw');
    pd(pd==0) = nan;
    
    candidateBlink_on = find(diff(isnan(pd)) == 1)+1;
    candidateBlink_off = find(diff(isnan(pd)) == -1)+1;
    
    L = min(length(candidateBlink_off), length(candidateBlink_on));
    candidateBlink_off = candidateBlink_off(1:L);
    candidateBlink_on = candidateBlink_on(1:L);
    candidateBlink_duration = (candidateBlink_off - candidateBlink_on)/samplingRate;
    
    blink_ind = (candidateBlink_duration > t_min) & (candidateBlink_duration < t_max);
    blink_on = candidateBlink_on(blink_ind);
    blink_off = candidateBlink_off(blink_ind);
    blink_duration = blink_off - blink_on;
    
    exptTime = range(gp.sub(sn).eye_data.tm)/60;
    
    blink = zeros(size(pd));
    blink(blink_off) = 1;
    gp.sub(sn).eye_data.blink = blink';
    blinkRate = length(blink_on) / exptTime;
    bR(sn) = blinkRate;
    fC(sn) = mean([gp.sub(sn).data.correct]);
    pD(sn) = nanmean(nanmean(gp.sub(sn).eye_data.pd));
end

%% autocorrelation of blinks
clear ac
for i = 1:120
    ac(i)  = corr(blink(1:end-i)', blink(1+i:end)');
end

figure(1); clf;
plot(ac)
ylim([-0.04 0.04])




%% compute average erp over all subjects for clicks on ====================
clear erp
erp{1} = compute_erp_basic_v1(gp, 'blink', 4, 0 , [-2 5], [-1 0], 'CLICKS_ON');
erp{2} = compute_erp_basic_v1(gp, 'pd', 4, 1 , [-2 5], [-1 0], 'CLICKS_ON');
erp{3} = compute_erp_basic_v1(gp, 'diff_pd', 4, 1 , [-2 5], [-1 0], 'CLICKS_ON');
%%
figure(1); clf; 
set(gcf, 'position', [0 0 500 600]);
ax = easy_gridOfEqualFigures([0.08 0.12 0.12 0.05], [0.13 0.03]);
axes(ax(1)); hold on;
[f,l,l0] = plot_erp_basic(gca, erp{1}, 30*60, 0);
plot([1 1], get(gca, 'ylim'),'k--')
set(f, 'facecolor', AZblue)
ylabel('blink rate [blinks/minute]')

axes(ax(2)); hold on;
[f,l] = plot_erp_basic(gca, erp{2}, 1, 0);
plot([1 1], get(gca, 'ylim'),'k--')
set(f, 'facecolor', AZred)
ylabel('pupil dilation [z-score]')

axes(ax(3)); hold on;
[f,l] = plot_erp_basic(gca, erp{3}, 1, 0);
set(f, 'facecolor', AZcactus)
plot([1 1], get(gca, 'ylim'),'k--')
ylabel('pupil asymmetry [z-score]')
set(ax, 'box', 'off')

% saveFigurePdf(gcf, '~/Desktop/clicks_basicErps')

%% erp based on agreement index ===========================================
clear erp l
erp{1} = compute_erp_ag_v1(gp, 'blink', 4, 0, [-10 2], [-1 0], 'CLICKS_ON');
erp{2} = compute_erp_ag_v1(gp, 'pd', 4, 1, [-10 2], [-1 0], 'CLICKS_ON');
erp{3} = compute_erp_ag_v1(gp, 'diff_pd', 4, 1, [-10 2], [-1 0], 'CLICKS_ON');
%%
clear M S
figure(1); clf; 
set(gcf, 'position', [0 0 500 600]);
ax = easy_gridOfEqualFigures([0.08 0.12 0.12 0.05], [0.13 0.03]);
axes(ax(1)); hold on;
[f,l, l0] = plot_erp_ag_v1(ax(1), erp{1}, 30*60, 0)
legend(f, {'ag=0' 'ag=1'}, 'location', 'northwest')
plot([1 1], get(gca, 'ylim'),'k--')
% set(f, 'facecolor', AZblue)
ylabel('blink rate [blinks/minute]')

axes(ax(2)); hold on;
[f,l, l0] = plot_erp_ag_v1(ax(2), erp{2}, 1, 0)
plot([1 1], get(gca, 'ylim'),'k--')
% set(f(1), 'facecolor', AZblue)
ylabel('pupil dilation [z-score]')

axes(ax(3)); hold on;
[f,l, l0] = plot_erp_ag_v1(ax(3), erp{3}, 1, 0)
plot([1 1], get(gca, 'ylim'),'k--')
% set(f, 'facecolor', AZblue)
ylabel('pupil asymmetry [z-score]')
% saveFigurePdf(gcf,'~/Desktop/clicks_agree')
%% look in t range for pupil
trange = [-10 -2];
for i = 1:2;%length(erp)
    Y = cat(3, erp{2}{i}.y_nb);
    tm = erp{1}{1}(1).tm';
    YY = squeeze(nanmean(Y,2));
    sig(i,:) = nanmean(YY( (tm>trange(1)) & (tm<trange(2)),:));
end

figure(1); clf; 
set(gcf, 'position', [0 0 500 300])
ax = easy_gridOfEqualFigures([0.15 0.1], [0.13 0.03]);
axes(ax); hold on;
plot(sig(1,:), sig(2,:),'o', 'markersize', 10, 'color', 'k', ...
    'markerfacecolor', AZsand)
plot([-0.4 0.3], [-0.4 0.3], 'k--')
[~,p] = ttest(sig(1,:), sig(2,:))
xlabel('pupil diameter for ag=0 [z-score]')
ylabel('pupil diameter for ag=1 [z-score]')
title(['average from ' num2str(trange(1)) 's to ' num2str(trange(2)) 's' ])
set(gca, 'tickdir', 'out')
saveFigurePdf(gcf, '~/Desktop/clicks_pupil_ag')
%% look in t range for pupil asymmetry
trange = [-10 0];
for i = 1:2;%length(erp)
    Y = cat(3, erp{3}{i}.y_nb);
    tm = erp{1}{1}(1).tm';
    YY = squeeze(nanmean(Y,2));
    sig(i,:) = nanmean(YY( (tm>trange(1)) & (tm<trange(2)),:));
end

figure(1); clf; 
set(gcf, 'position', [0 0 500 300])
ax = easy_gridOfEqualFigures([0.15 0.1], [0.13 0.03]);
axes(ax); hold on;
plot(sig(1,:), sig(2,:),'o', 'markersize', 10, 'color', 'k', ...
    'markerfacecolor', AZsand)
plot([-0.4 0.3], [-0.4 0.3], 'k--')
[~,p] = ttest(sig(1,:), sig(2,:))
xlabel('pupil asymmetry for ag=0 [z-score]')
ylabel('pupil asymmetry for ag=1 [z-score]')
title(['average from ' num2str(trange(1)) 's to ' num2str(trange(2)) 's' ])
set(gca, 'tickdir', 'out')
saveFigurePdf(gcf, '~/Desktop/clicks_pupil_asymmetry_ag2')


%% erp for error and correct ==============================================

clear erp l
% erp{1} = compute_erp_error_v1(gp, 'blink', 4, 0, [-10 2], [-1 0], 'CLICKS_ON');
% erp{2} = compute_erp_error_v1(gp, 'pd', 4, 1, [-10 2], [-1 0], 'CLICKS_ON');
% erp{3} = compute_erp_error_v1(gp, 'diff_pd', 4, 1, [-10 2], [-1 0], 'CLICKS_ON');
% erp{1} = compute_erp_error_v1(gp, 'blink', 4, 0, [-2 2], [-1 0], {'RESP_RIGHT' 'RESP_LEFT'});
% erp{2} = compute_erp_error_v1(gp, 'pd', 4, 1, [-2 2], [-1 0], {'RESP_RIGHT' 'RESP_LEFT'});
% erp{3} = compute_erp_error_v1(gp, 'diff_pd', 4, 1, [-2 2], [-1 0], {'RESP_RIGHT' 'RESP_LEFT'});
erp{1} = compute_erp_error_v1(gp, 'blink', 4, 0, [-2 2], [-1 0], {'OUT_CORRECT' 'OUT_WRONG'});
erp{2} = compute_erp_error_v1(gp, 'pd', 4, 1, [-2 2], [-1 0], {'OUT_CORRECT' 'OUT_WRONG'});
erp{3} = compute_erp_error_v1(gp, 'diff_pd', 4, 1, [-2 2], [-1 0], {'OUT_CORRECT' 'OUT_WRONG'});
%%
figure(1); clf; 
set(gcf, 'position', [0 0 500 600]);
ax = easy_gridOfEqualFigures([0.08 0.12 0.12 0.05], [0.13 0.03]);
axes(ax(1)); hold on;
[f,l, l0] = plot_erp_ag_v1(ax(1), erp{1}, 30*60, 0);
set(f(1), 'facecolor', [1 0.5 0.5])
set(f(2), 'facecolor', [0.5 0.5 1])
legend(f, {'error' 'correct'}, 'location', 'northwest')
plot([1 1], get(gca, 'ylim'),'k--')
% set(f, 'facecolor', AZblue)
ylabel('blink rate [blinks/minute]')

axes(ax(2)); hold on;
[f,l, l0] = plot_erp_ag_v1(ax(2), erp{2}, 1, 0);
set(f(1), 'facecolor', [1 0.5 0.5])
set(f(2), 'facecolor', [0.5 0.5 1])
plot([1 1], get(gca, 'ylim'),'k--')
% set(f(1), 'facecolor', AZblue)
ylabel('pupil dilation [z-score]')

axes(ax(3)); hold on;
[f,l, l0] = plot_erp_ag_v1(ax(3), erp{3}, 1, 0);
set(f(1), 'facecolor', [1 0.5 0.5])
set(f(2), 'facecolor', [0.5 0.5 1])

plot([1 1], get(gca, 'ylim'),'k--')
% set(f, 'facecolor', AZblue)
ylabel('pupil asymmetry [z-score]')
saveFigurePdf(gcf,'~/Desktop/clicks_error3')
%% look in t range for pupil
clear sig
trange = [-10 -2];
for i = 1:2;%length(erp)
    Y = cat(3, erp{2}{i}.y_nb);
    tm = erp{1}{1}(1).tm';
    YY = squeeze(nanmean(Y,2));
    sig(i,:) = nanmean(YY( (tm>trange(1)) & (tm<trange(2)),:));
end

figure(1); clf; 
set(gcf, 'position', [0 0 500 300])
ax = easy_gridOfEqualFigures([0.15 0.1], [0.13 0.03]);
axes(ax); hold on;
plot(sig(1,:), sig(2,:),'o', 'markersize', 10, 'color', 'k', ...
    'markerfacecolor', AZsand)
plot([-0.4 0.3], [-0.4 0.3], 'k--')
[~,p] = ttest(sig(1,:), sig(2,:))
xlabel('pupil diameter for error [z-score]')
ylabel('pupil diameter for correct [z-score]')
title(['average from ' num2str(trange(1)) 's to ' num2str(trange(2)) 's' ])
set(gca, 'tickdir', 'out')
saveFigurePdf(gcf, '~/Desktop/clicks_pupil_error')
%% look in t range for pupil asymmetry
clear sig
trange = [0 0.5];
for i = 1:2;%length(erp)
    Y = cat(3, erp{2}{i}.y_nb);
    tm = erp{1}{1}(1).tm';
    YY = squeeze(nanmean(Y,2));
    sig(i,:) = nanmean(YY( (tm>trange(1)) & (tm<trange(2)),:));
end

figure(1); clf; 
set(gcf, 'position', [0 0 500 300])
ax = easy_gridOfEqualFigures([0.15 0.1], [0.13 0.03]);
axes(ax); hold on;
plot(sig(1,:), sig(2,:),'o', 'markersize', 10, 'color', 'k', ...
    'markerfacecolor', AZsand)
plot([-1 0.3], [-1 0.3], 'k--')
ii = sig(1,:)>-1;
[~,p] = ttest(sig(1,ii), sig(2,ii))
xlabel('pupil asymmetry for error [z-score]')
ylabel('pupil asymmetry for correct [z-score]')
title(['average from ' num2str(trange(1)) 's to ' num2str(trange(2)) 's' ])
set(gca, 'tickdir', 'out')
saveFigurePdf(gcf, '~/Desktop/clicks_pupildiff_error')







%% average erp over subjects for out correct and out wrong ================
clear erp
for sn = 1:length(gp.sub)
    w = 4;
    tm = gp.sub(sn).eye_data.tm;
    pd = gp.sub(sn).eye_data.pd;
    for i = 1:size(pd,2)
        pd(:,i) = nanSmooth(pd(:,i), ones(w,1)/w, 'same');
    end
    pd = nanmean(pd,2);
    %pd = (diff(pd')');
    pd = (pd - nanmean(pd))/nanstd(pd);
    mess = gp.sub(sn).eye_mess;
    
    % clicks on
    sampling_rate = 30;
    T = [-10 2];
    BL = [-1 0];
    
    % ind = strcmp(mess.str, 'CLICKS_ON');
    % ind = strcmp(mess.str, 'RESP_LEFT') | strcmp(mess.str, 'RESP_RIGHT');
    ind = strcmp(mess.str, 'OUT_CORRECT'); 
    et = mess.tm(ind);
    iEvent = compute_iEvent(et, tm);
    erp1(sn) = compute_erp(iEvent, pd, T, BL, 'hi', sampling_rate);
    
    ind = strcmp(mess.str, 'OUT_WRONG');
    et = mess.tm(ind);
    iEvent = compute_iEvent(et, tm);
    erp0(sn) = compute_erp(iEvent, pd, T, BL, 'hi', sampling_rate);
    
end
%%
Y1 = cat(3, erp1.y_nb);
Y0 = cat(3, erp0.y_nb);
tm = erp1(1).tm';
YY1 = squeeze(nanmean(Y1,2));
M1 = nanmean(nanmean(Y1,3),2)*30*60;
S1 = nanstd(YY1, [], 2)/ sqrt(size(YY1, 2))*30*60;
YY0 = squeeze(nanmean(Y0,2));
M0 = nanmean(nanmean(Y0,3),2)*30*60;
S0 = nanstd(YY0, [], 2)/ sqrt(size(YY0, 2))*30*60;

figure(1); clf; hold on;
f = shadedErrorBars(tm', M1', S1');
l = plot(tm, M1, 'k', 'linewidth', 2);
set(f, 'facecolor', [0.5 0.5 1])
f(2) = shadedErrorBars(tm', M0', S0');
l = plot(tm, M0, 'k', 'linewidth', 2);
set(f(2), 'facecolor', [1 0.5 0.5])
yl = get(gca, 'ylim');
plot([0 0], yl, 'k--');
ylim(yl)
xlabel('time [seconds]')
ylabel('blink rate [blinks/minute]')
set(gca, 'tickdir', 'out')
legend(f, {'correct' 'wrong'})
% saveFigurePdf(gcf, '~/Desktop/blinks_rightWrong')
%%
trange = [0.5 1];
ind = (tm>trange(1)) & (tm<trange(2));
figure(1); clf; hold on;
x0 = nanmean(YY0(ind,:))*30*60;
x1 = nanmean(YY1(ind,:))*30*60;
ii = x0<inf;
plot(x0(ii), x1(ii),'.');
plot([0 120], [0 120],'k--')

[~,p] = ttest(x0(ii), x1(ii));

%% erp based on agreement index ===========================================
clear erp l
for sn = 1:length(gp.sub)
    data = gp.sub(sn).data;
    tm = gp.sub(sn).eye_data.tm;
    pd = gp.sub(sn).eye_data.pd;
    w = 10;
    for i = 1:size(pd,2)
        pd(:,i) = nanSmooth(pd(:,i), ones(w,1)/w, 'same');
    end
    pd = nanmean(pd,2);
    %pd = diff(pd')';
    pd = (pd - nanmean(pd))/nanstd(pd);
    mess = gp.sub(sn).eye_mess;
    
    sampling_rate = 30;
    T = [-20 2];
    BL = [-1 -0.5];
    
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
    
    
    ind = strcmp(mess.str, 'CLICKS_ON');
    %ind = strcmp(mess.str, 'RESP_RIGHT') | strcmp(mess.str, 'RESP_LEFT');
    et = mess.tm(ind);
    et = et(1:L);
    
    for i = 1:2
        ind = [data2.ag] == (i-1);
        iEvent = compute_iEvent(et(ind), tm);
        erp{i}(sn) = compute_erp(iEvent, pd, T, BL, 'hi', sampling_rate);
    end
    
    
end
%%
clear M S
for i = 1:2
    Y = cat(3, erp{i}.y_nb);
    tm = erp{1}(1).tm';
    YY = squeeze(nanmean(Y,2));
    M(:,i) = nanmean(nanmean(Y,3),2);
    S(:,i) = nanstd(YY, [], 2)/ sqrt(size(YY, 2));
    
end
figure(1); clf; hold on;
for i = 1:size(M,2)
    f(i) = shadedErrorBars(tm', M(:,i)', S(:,i)');
    l(i) = plot(tm, M(:,i), 'k', 'linewidth', 2);
end
set(f(1), 'facecolor', [1 0.5 0])
set(f(2), 'facecolor', [0 0.7 1])

yl = get(gca, 'ylim');
plot([0 0], yl, 'k--');
ylim(yl)
xlabel('time [seconds]')
ylabel('pupil diameter')
set(gca, 'tickdir', 'out')
%% look in t range
trange = [-10 -2];
for i = 1:length(erp)
    Y = cat(3, erp{i}.y_nb);
    tm = erp{1}(1).tm';
    YY = squeeze(nanmean(Y,2));
    sig(i,:) = nanmean(YY( (tm>trange(1)) & (tm<trange(2)),:));
end

figure(1); clf; hold on;
plot(sig(1,:), sig(2,:),'.')
plot([-0.4 0.3], [-0.4 0.3])
[~,p] = ttest(sig(1,:), sig(2,:))

%% erp based on mistake








%%
clear l
figure(1); clf; 
ax = easy_gridOfEqualFigures([0.1 0.1 0.03], [0.1 0.03]);
axes(ax(1)); hold on;
for i = 1:length(erp)
    l(i) = plot(erp(i).tm, (erp(i).y'));
end

set(l(1), 'color', 'r')
set(l(2), 'color', 'g')




%% compute erp
sn = 20;
w = 4;
tm = gp.sub(sn).eye_data.tm;
pd = gp.sub(sn).eye_data.blink;
for i = 1:size(pd,2)
    pd(:,i) = nanSmooth(pd(:,i), ones(w,1)/w, 'same');
end
mess = gp.sub(sn).eye_mess;

% clicks on
sampling_rate = 30;
T = [-2 2];
BL = [-1 0];

ind = strcmp(mess.str, 'CLICKS_ON');
% ind = strcmp(mess.str, 'RESP_LEFT') | strcmp(mess.str, 'RESP_RIGHT');
% ind = strcmp(mess.str, 'OUT_CORRECT') | strcmp(mess.str, 'OUT_WRONG');
et = mess.tm(ind);
iEvent = compute_iEvent(et, tm);
erp = compute_erp(iEvent, pd, T, BL, 'hi', sampling_rate);

figure(1); clf;
% ax = easy_gridOfEqualFigures([0.1 0.1 0.05], [0.1 0.03]);
% axes(ax(1)); hold on;
plot(erp.tm, (erp.y_nb'))

% axes(ax(2)); hold on;
% plot(erp.tm, diff(erp.y'))

%% response left vs right
sn = 11;
w = 10;
tm = gp.sub(sn).eye_data.tm;
pd = gp.sub(sn).eye_data.pd;
for i = 1:size(pd,2)
    pd(:,i) = nanSmooth(pd(:,i), ones(w,1)/w, 'same');
end
mess = gp.sub(sn).eye_mess;

T = [-2 2];
BL = [-1 0];

ind = strcmp(mess.str, 'OUT_CORRECT');
et = mess.tm(ind);
iEvent = compute_iEvent(et, tm);
erpL = compute_erp(iEvent, pd, T, BL, 'hi', sampling_rate);

ind = strcmp(mess.str, 'OUT_WRONG');
et = mess.tm(ind);
iEvent = compute_iEvent(et, tm);
erpR = compute_erp(iEvent, pd, T, BL, 'hi', sampling_rate);

figure(1); clf;
% ax = easy_gridOfEqualFigures([0.1 0.1 0.05], [0.1 0.03]);
% axes(ax(1)); 
hold on;
plot(erpL.tm, mean(erpL.y_nb,2))
plot(erpR.tm, mean(erpR.y_nb,2))

% axes(ax(2)); hold on;
% plot(erpL.tm, diff(erpL.y'))
% plot(erpR.tm, diff(erpR.y'))

%% correct vs wrong
sn = 11;
w = 5;
tm = gp.sub(sn).eye_data.tm;
pd = gp.sub(sn).eye_data.pd;
for i = 1:size(pd,2)
    pd(:,i) = nanSmooth(pd(:,i), ones(w,1)/w, 'same');
end
mess = gp.sub(sn).eye_mess;

data = gp.sub(sn).data;

T = [-1 2];
BL = [-1 -0];

DN = [2 4 6 8 10 12];
dn = abs([data.nL] - [data.nR]);
ind = strcmp(mess.str, 'CLICKS_ON');
% ind = strcmp(mess.str, 'RESP_RIGHT') | strcmp(mess.str, 'RESP_LEFT');

et = mess.tm(ind);
et = et(([data(1:length(et)).correct] == 1) & ismember(dn(1:length(et)), DN));
iEvent = compute_iEvent(et, tm);
erp1 = compute_erp(iEvent, pd, T, BL, 'hi', sampling_rate);

% ind = strcmp(mess.str, 'RESP_RIGHT') | strcmp(mess.str, 'RESP_LEFT');
ind = strcmp(mess.str, 'CLICKS_ON');
et = mess.tm(ind);
et = et(([data(1:length(et)).correct] == 0)  & ismember(dn(1:length(et)), DN));
iEvent = compute_iEvent(et, tm);
erp0 = compute_erp(iEvent, pd, T, BL, 'hi', sampling_rate);

figure(1); clf;
% ax = easy_gridOfEqualFigures([0.1 0.1 0.05], [0.1 0.03]);
% axes(ax(1)); 
hold on;
plot(erp1.tm, mean(erp1.y_nb,2), 'r')
plot(erp0.tm, mean(erp0.y_nb,2), 'b')

% axes(ax(2)); hold on;
% plot(erp1.tm, diff(erp1.y'), 'r')
% plot(erp0.tm, diff(erp0.y'), 'b')

%% clicks on by abs(nL-nR)
sn = 4;
data = gp.sub(sn).data;
tm = gp.sub(sn).eye_data.tm;
pd = gp.sub(sn).eye_data.pd;
for i = 1:size(pd,2)
    pd(:,i) = nanSmooth(pd(:,i), ones(1,1)/1, 'same');
end
mess = gp.sub(sn).eye_mess;

clear erp l1 l2 l

T = [-4 2];
BL = [-1 -0.5];



% ind = strcmp(mess.str, 'CLICKS_ON');
% ind = strcmp(mess.str, 'RESP_RIGHT') | strcmp(mess.str, 'RESP_LEFT');
ind = strcmp(mess.str, 'OUT_CORRECT') | strcmp(mess.str, 'OUT_WRONG');
et = mess.tm(ind);

dn = abs([data.nL] - [data.nR]);
df = [0:2:6];
% df(1) = [0 2 8 ];
% df(2) = [0  

for i = 1:length(df)
    ind = (dn==df(i));
    iEvent = compute_iEvent(et(ind), tm);
    erp(i) = compute_erp(iEvent, pd, T, BL, 'hi', sampling_rate);
end


figure(1); clf;
ax = easy_gridOfEqualFigures([0.1 0.1 0.05], [0.1 0.03]);
axes(ax(1)); hold on;
for i = 1:length(erp)
    l1(i) = plot(erp(i).tm, mean(erp(i).y'));
end

axes(ax(2)); hold on;
for i = 1:length(erp)
    l2(i) = plot(erp(i).tm, diff(erp(i).y'));
end
for i = 1:length(erp)
    col{i} = colormix([1 0 0], [0 0 1], i, length(erp));
    set([l1(i) l2(i)], 'color', col{i});
end

















%% choice triggered covariance
i = 1;
% left = 1
goL = double([data([data.choice]==1).lClicks]);
goR = double([data([data.choice]==2).lClicks]);
















%% psychometric functions =================================================
%% basic choice curve and RT curve
clear m s mRT sRT
for sn = 1:length(gp.sub)
    d_vals = [-14:2:14];
    [m(sn,:), s(sn,:), mRT(sn,:), sRT(sn,:), mAG(sn,:)] = compute_psychometricCurve(gp.sub(sn).data, d_vals, 5);
end
M = nanmean(m);
S = nanstd(m) / sqrt(size(m,1));
MRT = nanmean(mRT);
SRT = nanstd(mRT) / sqrt(size(mRT,1));
MAG = nanmean(mAG);
SAG = nanstd(mAG) / sqrt(size(mAG,1));

%% compute expected agreement index given m
tAG = m.*m + (1-m).*(1-m);
TAG = nanmean(tAG);
StAG = nanstd(tAG) / sqrt(size(tAG,1));

%%
figure(1); clf; 
set(gcf, 'position', [0 0 500 600]);
ax = easy_gridOfEqualFigures([0.12 0.13 0.13 0.05], [0.13 0.03]);

axes(ax(1)); hold on;
plot(d_vals, m, 'color', [1 1 1]*0.75)
errorbar(d_vals, M, S, 'color', 'k', 'linewidth', 5)
xlabel('n_{right} - n_{left}')
ylabel('p_{right}')

axes(ax(2)); hold on;
plot(d_vals, mRT, 'color', [1 1 1]*0.75)
errorbar(d_vals, MRT, SRT, 'color', 'k', 'linewidth', 5)
xlabel('n_{right} - n_{left}')
ylabel('RT [s]')
set(ax, 'box', 'off', 'tickdir', 'out')

axes(ax(3)); hold on;
plot(d_vals, mAG, 'color', [1 1 1]*0.75)
errorbar(d_vals, MAG, SAG, 'color', 'k', 'linewidth', 5)
errorbar(d_vals, TAG, StAG, 'color', 'r', 'linewidth', 3)
xlabel('n_{right} - n_{left}')
ylabel('p_{agree}')
set(ax, 'box', 'off', 'tickdir', 'out')
saveFigurePdf(gcf, '~/Desktop/clicks_behavior_ag')


 
    
    
    
    
    
    
%% median split on RT
clear m s mRT sRT
for sn = 1:length(gp.sub)
    pd = gp.sub(sn).eye_data.pd;
    pd = gp.sub(sn).eye_mess
    
    [m(:,:,sn), s, mRT(:,:,sn), sRT] = compute_psychometricCurve_medianSplitRTs(...
        gp.sub(sn).data, d_vals, 5)
end
M = nanmean(m,3);
S = nanstd(m,[],3) / sqrt(size(m,3));
MRT = nanmean(mRT,3);
SRT = nanstd(mRT,[],3) / sqrt(size(mRT,3));

%%
figure(1); clf
ax = easy_gridOfEqualFigures([0.1 0.15 0.05], [0.15 0.03]);
axes(ax(1)); hold on;
errorbar([d_vals' d_vals'], M, S)
xlabel('n_{right} - n_{left}')
ylabel('p_{right}')
axes(ax(2));
errorbar([d_vals' d_vals'], MRT, SRT)
xlabel('n_{right} - n_{left}')
ylabel('RT [s]')
set(ax, 'box', 'off', 'tickdir', 'out')
saveFigurePdf(gcf, '~/Desktop/clicks_behavior_medSplitRT')


%% median split on pupil
clear erp
erp = compute_erp_basic_v1(gp, 'pd', 1, 1 , [-20 2], [-1 0], 'CLICKS_ON');
trange = [-10 -2];
%% 
clear m mRT mAG
for sn = 1:length(erp)
    ind = (erp(sn).tm > trange(1)) & (erp(sn).tm < trange(2));
    sig = squeeze(nanmean(erp(sn).D(ind,:,:),1))';
    % focus on first 500 trials
    L = min(length(sig), 500);
    sig = sig(1:L);
    
    [m(:,:,sn), s, mRT(:,:,sn), sRT, mAG(:,:,sn)] = compute_psychometricCurve_medianSplitOnX_v1(...
        gp.sub(sn).data(1:L), d_vals, sig, 10);
end
    
M = nanmean(m,3);
S = nanstd(m,[],3) / sqrt(size(m,3));
MRT = nanmean(mRT,3);
SRT = nanstd(mRT,[],3) / sqrt(size(mRT,3));
MAG = nanmean(mAG,3);
SAG = nanstd(mAG,[],3) / sqrt(size(mAG,3));

figure(1); clf
ax = easy_gridOfEqualFigures([0.1 0.15 0.15 0.05], [0.15 0.03]);
axes(ax(1)); hold on;
errorbar([d_vals' d_vals'], M, S)
xlabel('n_{right} - n_{left}')
ylabel('p_{right}')
legend({'low pupil' 'high pupil'}, 'location', 'northwest')
axes(ax(2));
errorbar([d_vals' d_vals'], MRT, SRT)
xlabel('n_{right} - n_{left}')
ylabel('RT [s]')
set(ax, 'box', 'off', 'tickdir', 'out')
axes(ax(3));
errorbar([d_vals' d_vals'], MAG, SAG)
xlabel('n_{right} - n_{left}')
ylabel('p_{agree}')
set(ax, 'box', 'off', 'tickdir', 'out')
saveFigurePdf(gcf, '~/Desktop/clicks_behavior_medSplitPupil_ag')
   

%% scatter plot of agreement split by pupil
xx = squeeze(nanmean(mAG,1));
figure(1); clf; hold on;
plot(xx(1,:), xx(2,:),'o', ...
    'markersize', 10, ...
    'markerfacecolor', AZsand, ...
    'color', 'k');
plot([0 1], [0 1], 'k--')
[~,p] = ttest(xx(1,:), xx(2,:));
xlabel('average agreement index, low pupil')
ylabel('average agreement index, high pupil')
% saveFigurePdf(gcf, '~/Desktop/pupilsplit_agreement')

%% scatter plot of agreement split by pupil
xx = squeeze(nanmean(mRT,1));
figure(1); clf; hold on;
plot(xx(1,:), xx(2,:),'o', ...
    'markersize', 10, ...
    'markerfacecolor', AZsand, ...
    'color', 'k');
plot([1 2.5],[1 2.5], 'k--')
[~,p] = ttest(xx(1,:), xx(2,:));
xlabel('reaction time, low pupil')
ylabel('reaction time, high pupil')
saveFigurePdf(gcf, '~/Desktop/pupilsplit_RT')

%% p low
p_low = squeeze(nanmean(m(d_vals'<0,:,:)+(1-m(d_vals'>0,:,:)),1))/2;
figure(1); clf; hold on;

plot(p_low(1,:), p_low(2,:),'o', ...
    'markersize', 10, ...
    'markerfacecolor', AZsand, ...
    'color', 'k');
plot([0 1],[0 1], 'k--')
[~,p] = ttest(p_low(1,:), p_low(2,:));
xlabel('p_{low}, low pupil')
ylabel('p_{low}, high pupil')
saveFigurePdf(gcf, '~/Desktop/plow')






%% choice-triggered average of clicks =====================================
for sn = 1:length(gp.sub)
    data = gp.sub(sn).data;
    L = min(length(data), 500);
    for i = 1:2
        ind = find([data(1:L).choice] == i);
        lClicks(:,i,sn) = mean([data(ind).lClicks],2);
        rClicks(:,i,sn) = mean([data(ind).rClicks],2);
    end
end

lM = nanmean(lClicks,3);
lS = nanstd(lClicks,[],3)/sqrt(size(lClicks,3));
figure(1); clf; hold on;
l = errorbar(lM, lS);
ylim([0 1]);
set(l(1), 'color', AZblue, 'linewidth', 2)
set(l(2), 'color', AZred, 'linewidth', 2)
xlim([0 21])
xlabel('time step')
ylabel('probability of left click')
legend({'left choice' 'right choice'})
% saveFigurePdf(gcf, '~/Desktop/clickTriggeredAverage')


%% choice-triggered average of clicks split out differnet dns
for sn = 1:length(gp.sub)
    data = gp.sub(sn).data;
    L = min(length(data), 500);
    dn = [data.nR]-[data.nL];
    for i = 1:2
        
        ind = find(([data(1:L).choice] == i)&(abs(dn(1:L))==0));
        lClicks(:,i,sn) = mean([data(ind).lClicks],2);
        rClicks(:,i,sn) = mean([data(ind).rClicks],2);
    end
end

lM = nanmean(lClicks,3);
lS = nanstd(lClicks,[],3)/sqrt(size(lClicks,3));
figure(1); clf; hold on;
l = errorbar(lM, lS);
ylim([0 1]);
set(l(1), 'color', AZblue, 'linewidth', 2)
set(l(2), 'color', AZred, 'linewidth', 2)
xlim([0 21])
xlabel('time step')
ylabel('probability of left click')
legend({'left choice' 'right choice'})
% saveFigurePdf(gcf, '~/Desktop/clickTriggeredAverage_dn')

%% choice-triggered average of clicks median split on RT
for sn = 1:length(gp.sub)
    data = gp.sub(sn).data;
    L = min(length(data), 500);
    dn = [data.nR]-[data.nL];
    RT = [data.RT];
    RT = RT(1:L);
    for i = 1:2
        
        
        ind = (([data(1:L).choice] == i));%&(dn(1:L)==4));
        rt = RT(ind);
        lClicks(:,i,sn) = mean([data(ind & RT<median(rt)).lClicks],2);
        rClicks(:,i,sn) = mean([data(ind).rClicks],2);
    end
end

lM = nanmean(lClicks,3);
lS = nanstd(lClicks,[],3)/sqrt(size(lClicks,3));
dM = -nanmean(lClicks(:,2,:)-lClicks(:,1,:),3);
dS = nanstd(lClicks(:,2,:)-lClicks(:,1,:),[],3)/sqrt(size(lClicks,3));
figure(1); clf; 
ax = easy_gridOfEqualFigures([0.1 0.1 0.1], [0.1 0.03]);
axes(ax(1)); hold on;
l = errorbar(lM, lS);
ylim([0 1]);
set(l(1), 'color', AZblue, 'linewidth', 2)
set(l(2), 'color', AZred, 'linewidth', 2)
xlim([0 21])
xlabel('time step')
ylabel('probability of left click')
legend({'left choice' 'right choice'})
% saveFigurePdf(gcf, '~/Desktop/clickTriggeredAverage_dn')

axes(ax(2)); hold on;
l = errorbar(dM, dS);
ylim([0 0.3]);
set(l(1), 'color', AZblue, 'linewidth', 2)
xlim([0 21])
xlabel('time step')
ylabel('probability of left click')
% saveFigurePdf(gcf, '~/Desktop/clickTriggeredAverage_dn')









