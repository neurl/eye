function zz0 = nanSmooth(xx, aa, opt)

% pd_filtered = nanSmooth(pd, ones(60,1), 'same');

nn = conv(double(~isnan(xx)), aa, opt);

xx0 = xx;
xx0(isnan(xx)) = 0;

yy0 = conv(xx0, aa, opt);
zz0 = yy0 ./ nn;
