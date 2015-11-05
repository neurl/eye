global orange fontname fontsize


% fontname = 'arial';
% fontsize = 12;
% ABCfontsize = 20;
% fontweight = 'bold';
% linewidth = 2;

fontname = 'arial';
fontsize = 12;
ABCfontsize = 24;
fontweight = 'normal';
linewidth = 2;


orange = [0.906 0.463 0.247];


% colormap gray
% CC = colormap;
% CM = (CC).*repmat((1-[0.906 0.463 0.247]/0.906), size(CC,1),1);
% colormap(1-CM)

set(0, 'defaultfigurecolor', 'w', ...
    'defaultaxesfontsize', fontsize, ...
    'defaultaxesfontweight', fontweight, ...
    'defaultaxesfontname', fontname, ...
    'defaultaxestickdir', 'out', ...
    'defaultaxesbox', 'on', ...
    'defaultaxesydir', 'normal', ...
    'defaultlinelinewidth', linewidth, ...
    'defaultlinemarkersize', 20)