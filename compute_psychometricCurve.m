function [m, s, mRT, sRT] = compute_psychometricCurve(data, d_vals)

df = [data.nR] - [data.nL];

for i = 1:length(d_vals)
    ind = df == d_vals(i);
    
    % choice
    X = [data(ind).choice]==2;
    m(i) = mean(X);
    alpha = sum(X);
    beta = sum(~X);
    s(i) = sqrt(alpha * beta / (alpha + beta)^2 / (alpha + beta + 1));
    
    % RT
    Y = [data(ind).RT];
    mRT(i) = mean(Y);
    sRT(i) = std(Y) / sqrt(sum(ind));
end
