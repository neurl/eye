function [f,l, l0] = plot_erp_ag_v1(ax, erp, sf, blflag)


for i = 1:2
    switch blflag
        case 1
            Y = cat(3, erp{i}.y)*sf;
        case 0
            Y = cat(3, erp{i}.y_nb)*sf;
    end
    tm = erp{1}(1).tm';
    YY = squeeze(nanmean(Y,2));
    M(:,i) = nanmean(nanmean(Y,3),2);
    S(:,i) = nanstd(YY, [], 2)/ sqrt(size(YY, 2));
    
end
axes(ax); hold on;
for i = 1:size(M,2)
    f(i) = shadedErrorBars(tm', M(:,i)', S(:,i)');
    l(i) = plot(tm, M(:,i), 'k', 'linewidth', 2);
end
set(f(1), 'facecolor', [1 0.5 0])
set(f(2), 'facecolor', [0 0.7 1])

yl = get(ax, 'ylim');
l0 = plot([0 0], yl, 'k--');
ylim(yl)
xlabel('time [seconds]')
ylabel('pupil diameter')
set(ax, 'tickdir', 'out')