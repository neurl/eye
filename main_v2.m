clear

maindir = ['~/Work/Projects/clicks/task_clicks/'];
datadir = [maindir 'datadir/2015_frozenNoise_v2/'];
fundir  = [maindir 'fundir/v1/'];

addpath(fundir);
cd(fundir);

gp = group(datadir);
gp.load;


%% clean pupil
gp.clean_pupil;


%% detect blinks
for sn = 1:length(gp.sub)
    t_min = 50/1000; t_max = 500/1000;
    
    samplingRate = 1/(gp.sub(1).eye_data.tm(2)-gp.sub(1).eye_data.tm(1));
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
    gp.sub(sn).eye_data.blink = blink;
    blinkRate = length(blink_on) / exptTime;
    bR(sn) = blinkRate;
    fC(sn) = mean([gp.sub(sn).data.correct]);
    pD(sn) = nanmean(nanmean(gp.sub(sn).eye_data.pd));
end

%% autocorrelation of blinks
clear ac
for i = 1:120
    ac(i)  = corr(pd(1:end-i), pd(1+i:end));
end

figure(1); clf;
plot(ac)
ylim([-0.04 0.04])


%% compute erp
sn = 4;
w = 4;
tm = gp.sub(sn).eye_data.tm;
pd = gp.sub(sn).eye_data.blink';
for i = 1:size(pd,2)
    pd(:,i) = nanSmooth(pd(:,i), ones(w,1)/w, 'same');
end
mess = gp.sub(sn).eye_mess;

% clicks on
sampling_rate = 60;
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
pd = gp.sub(sn).eye_data.blink';
for i = 1:size(pd,2)
    pd(:,i) = nanSmooth(pd(:,i), ones(w,1)/w, 'same');
end
mess = gp.sub(sn).eye_mess;

sampling_rate = 60;
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
pd = gp.sub(sn).eye_data.blink';
for i = 1:size(pd,2)
    pd(:,i) = nanSmooth(pd(:,i), ones(w,1)/w, 'same');
end
mess = gp.sub(sn).eye_mess;

data = gp.sub(sn).data;
sampling_rate = 60;
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
sampling_rate = 60;
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

%% erp based on agreement index
sn = 7;
data = gp.sub(sn).data;
tm = gp.sub(sn).eye_data.tm;
pd = gp.sub(sn).eye_data.pd;
for i = 1:size(pd,2)
    pd(:,i) = nanSmooth(pd(:,i), ones(1,1)/1, 'same');
end
mess = gp.sub(sn).eye_mess;

clear erp
sampling_rate = 60;
T = [-1 2];
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
% ind = strcmp(mess.str, 'RESP_RIGHT') | strcmp(mess.str, 'RESP_LEFT');
et = mess.tm(ind);
et = et(1:L);



for i = 1:2
    ind = [data2.ag] == (i-1);
    iEvent = compute_iEvent(et(ind), tm);
    erp(i) = compute_erp(iEvent, pd, T, BL, 'hi', sampling_rate);
end
clear l
figure(1); clf; 
ax = easy_gridOfEqualFigures([0.1 0.1 0.03], [0.1 0.03]);
axes(ax(1)); hold on;
for i = 1:length(erp)
    l(i) = plot(erp(i).tm, mean(erp(i).y'));
end

set(l(1), 'color', 'r')
set(l(2), 'color', 'g')

axes(ax(2)); hold on;
for i = 1:length(erp)
    l(i) = plot(erp(i).tm, diff(erp(i).y'));
end

set(l(1), 'color', 'r')
set(l(2), 'color', 'g')







%% choice-triggered average of clicks
for i = 1:2
    ind = [data.choice] == i;
    
    lClicks(:,i) = mean([data(ind).lClicks],2);
    rClicks(:,i) = mean([data(ind).rClicks],2);
end




%% choice triggered covariance
i = 1;
% left = 1
goL = double([data([data.choice]==1).lClicks]);
goR = double([data([data.choice]==2).lClicks]);
















%% psychometric functions =================================================
%% basic choice curve and RT curve
sn = 14;
d_vals = [-14:2:14];
[m, s, mRT, sRT] = compute_psychometricCurve(gp.sub(sn).data, d_vals);

figure(1); clf; 
ax = easy_gridOfEqualFigures([0.1 0.15 0.05], [0.15 0.03]);
axes(ax(1)); hold on;
errorbar(d_vals, m, s)
xlabel('n_{right} - n_{left}')
ylabel('p_{right}')
axes(ax(2));
errorbar(d_vals, mRT, sRT)
xlabel('n_{right} - n_{left}')
ylabel('RT [s]')
set(ax, 'box', 'off', 'tickdir', 'out')

%% median split on RT
sn = 7;
[m, s, mRT, sRT] = compute_psychometricCurve_medianSplitRTs(gp.sub(sn).data, d_vals);
figure(1); clf
ax = easy_gridOfEqualFigures([0.1 0.15 0.05], [0.15 0.03]);
axes(ax(1)); hold on;
errorbar([d_vals' d_vals'], m, s)
xlabel('n_{right} - n_{left}')
ylabel('p_{right}')
axes(ax(2));
errorbar([d_vals' d_vals'], mRT, sRT)
xlabel('n_{right} - n_{left}')
ylabel('RT [s]')
set(ax, 'box', 'off', 'tickdir', 'out')






