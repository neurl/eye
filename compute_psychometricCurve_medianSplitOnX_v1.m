function [m, s, mRT, sRT, mAG, sAG] = compute_psychometricCurve_medianSplitOnX_v1(data, d_vals, XX, RTmax)

% focus on first 500 (if we have 500)
L = min(500, length(data));
L = L - mod(L,2);
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
    
    D = data2(ind);
    
    RT = [data2(ind).RT];
    X = XX(ind);
    i1 = X<median(X);
    i2 = X>=median(X);
    
    % choice
    X1 = [D(i1).choice]==2;
    X2 = [D(i2).choice]==2;
    
    m(i,1) = mean(X1);
    m(i,2) = mean(X2);
    alpha = sum(X1);
    beta = sum(~X1);
    s(i,1) = sqrt(alpha * beta / (alpha + beta)^2 / (alpha + beta + 1));
    
    alpha = sum(X2);
    beta = sum(~X2);
    s(i,2) = sqrt(alpha * beta / (alpha + beta)^2 / (alpha + beta + 1));
    
    % RT
    Y1 = [D(i1).RT];
    Y1 = Y1(Y1<RTmax);
    mRT(i,1) = mean(Y1);
    sRT(i,1) = std(Y1) / sqrt(sum(i1));
    
    Y2 = [D(i2).RT];
    Y2 = Y2(Y2<RTmax);
    mRT(i,2) = mean(Y2);
    sRT(i,2) = std(Y2) / sqrt(sum(i2));
    
    % agreement
    Z1 = [D(i1).ag];
    mAG(i,1) = mean(Z1);
    sAG(i,1) = std(Z1) / sqrt(sum(ind));
    
    Z2 = [D(i2).ag];
    mAG(i,2) = mean(Z2);
    sAG(i,2) = std(Z2) / sqrt(sum(ind));
    
end
