clear

maindir = ['~/Dropbox/NeuRLLab_git/eye-data/'];
datadir = [maindir 'datadir/'];
fundir  = [maindir 'fundir/'];

addpath(fundir);
cd(fundir);

cd(datadir);
d = dir('*.mat');

sn = 1;
beh_name = d(sn).name;
eye_name = [d(sn).name(1:end-4) '_eye.tsv'];
load(beh_name);

% data = data(1:500);



%% load eye data
eye_data_files = dir([datadir '*_eye.tsv']);
importeyedata() 

% fid = fopen(eye_name);
% 
% X = textscan(fid, '%s', 'delimiter', '\n');
% fclose(fid);
% T = textscan(X{1}{1}, '%s', 'delimiter', '\t');
% 
% 
% Y = strvcat(X{1}{:});
% ind_mess = Y(:,1) == 'M';
% M = {X{1}{ind_mess}};
% Z = {X{1}{~ind_mess}};
% Z = {Z{2:end}};
% A = strvcat(Z{:});
% B = textscan(A','%s%s%f%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f');
% 
% E = strvcat(B{2});
% sc = str2num(E(:,end-5:end));
% mn = str2num(E(:,end-8:end-7))*60;
% hr = str2num(E(:,end-11:end-10))*60*60;
% tm = hr+mn+sc;
% 
% eye_data.tm = tm;
% eye_data.pd_raw = [B{15} B{22}];
% 
% eye_data.var_names = {T{1}{4:end}};
% 
% 
% % messages
% M = {M{2:end}};
% C = strvcat(M{:});
%     
% D = textscan(C','%s%s%s%s%s');
% F = strvcat(D{3});
% sc = str2num(F(:,end-5:end));
% mn = str2num(F(:,end-8:end-7))*60;
% hr = str2num(F(:,end-11:end-10))*60*60;
% tm = hr+mn+sc;
% 
% 
% mess.tm = tm;
% mess.str = D{5};


%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s

% % choice-triggered average of clicks
% for i = 1:2
%     ind = [data.choice] == i;
%     
%     lClicks(:,i) = mean([data(ind).lClicks],2);
%     rClicks(:,i) = mean([data(ind).rClicks],2);
% end
% 
% % psychometric functions =================================================
% % basic choice curve and RT curve
% d_vals = [-14:2:14];
% [m, s, mRT, sRT] = compute_psychometricCurve(data, d_vals);
% 
% figure(1); clf; 
% ax = easy_gridOfEqualFigures([0.1 0.15 0.05], [0.15 0.03]);
% axes(ax(1)); hold on;
% errorbar(d_vals, m, s)
% xlabel('n_{right} - n_{left}')
% ylabel('p_{right}')
% axes(ax(2));
% errorbar(d_vals, mRT, sRT)
% xlabel('n_{right} - n_{left}')
% ylabel('RT [s]')
% set(ax, 'box', 'off', 'tickdir', 'out')
% 
% % median split on RT
% [m, s, mRT, sRT] = compute_psychometricCurve_medianSplitRTs(data, d_vals);
% figure(1); clf
% ax = easy_gridOfEqualFigures([0.1 0.15 0.05], [0.15 0.03]);
% axes(ax(1)); hold on;
% errorbar([d_vals' d_vals'], m, s)
% xlabel('n_{right} - n_{left}')
% ylabel('p_{right}')
% axes(ax(2));
% errorbar([d_vals' d_vals'], mRT, sRT)
% xlabel('n_{right} - n_{left}')
% ylabel('RT [s]')
% set(ax, 'box', 'off', 'tickdir', 'out')
% 
% 
% 
% 
% 
% 
