function [f,l, l0] = plot_erp_basic(ax, erp, sf, blflag)

switch blflag
    case 1
        Y = cat(3, erp.y)*sf;
    case 0
        Y = cat(3, erp.y_nb)*sf;
end
tm = erp(1).tm';
YY = squeeze(nanmean(Y,2));
M = nanmean(nanmean(Y,3),2);
S = nanstd(YY, [], 2)/ sqrt(size(YY, 2));

axes(ax); hold on;
f = shadedErrorBars(tm', M', S');
l = plot(tm, M, 'k', 'linewidth', 2);
set(f, 'facecolor', [1 0.7 1])
yl = get(gca, 'ylim');
l0 = plot([0 0], yl, 'k--');
ylim(yl)
xlabel('time [seconds]')
ylabel('pupil diameter')
set(gca, 'tickdir', 'out')
