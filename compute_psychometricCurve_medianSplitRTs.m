function [m, s, mRT, sRT] = compute_psychometricCurve_medianSplitRTs(data, d_vals)

df = [data.nR] - [data.nL];

for i = 1:length(d_vals)
    ind = df == d_vals(i);
    
    D = data(ind);
    
    RT = [data(ind).RT];
    i1 = RT<median(RT);
    i2 = RT>=median(RT);
    
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
    mRT(i,1) = mean(Y1);
    sRT(i,1) = std(Y1) / sqrt(sum(i1));
    
    Y2 = [D(i2).RT];
    mRT(i,2) = mean(Y2);
    sRT(i,2) = std(Y2) / sqrt(sum(i2));
end
