function [m, s, mRT, sRT, mAG, sAG] = compute_psychometricCurve(data, d_vals, RTmax)

% also compute agreement index



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

df = [data2.nR] - [data2.nL];


for i = 1:length(d_vals)
    ind = df == d_vals(i);
    
    % choice
    X = [data2(ind).choice]==2;
    m(i) = mean(X);
    alpha = sum(X);
    beta = sum(~X);
    s(i) = sqrt(alpha * beta / (alpha + beta)^2 / (alpha + beta + 1));
    
    % RT
    Y = [data2(ind).RT];
    Y = Y(Y<RTmax);
    mRT(i) = mean(Y);
    sRT(i) = std(Y) / sqrt(sum(ind));
    
    % agreement
    Y = [data2(ind).ag];
    mAG(i) = mean(Y);
    sAG(i) = std(Y) / sqrt(sum(ind));
    
end
