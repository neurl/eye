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
global AZred AZblue AZcactus AZsky AZriver AZsand AZmesa AZbrick

AZred = [171,5,32]/256;
AZblue = [12,35,75]/256;
AZcactus = [92, 135, 39]/256;
AZsky = [132, 210, 226]/256;
AZriver = [7, 104, 115]/256;
AZsand = [241, 158, 31]/256;
AZmesa = [183, 85, 39]/256;
AZbrick = [74, 48, 39]/256;

% lighten blue
AZblue = 1*AZblue + 0*[1 1 1];

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