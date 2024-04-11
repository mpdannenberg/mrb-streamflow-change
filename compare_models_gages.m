
T = readtable('./data/Missouri_Big6_UpdatedUnregulatedFlow_DailyMonthlyWY_avgCMS.xlsx',...
    'Sheet','WaterYear_cms');
T.Properties.VariableNames = {'Year','FTRA','SUX','NCNE','MKC','MNMO','HEMO'};

basins = shaperead('./data/LMRBgagebasins/lmrb_shapefiles.shp');

%% Make figure showing comparison between gage data and MsTMIP
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6.5 3.];

load ./data/MsTMIP_WaterBudget;
[~, ia, ib] = intersect(year, T.Year);

% HEMO - unscaled
Q = WB_SG3_unscaled;

subplot(2,4,[1 3])
plot(year(ia)', Q(ia,:), '-', 'Color',[0.8 0.8 0.8])
hold on;
plot(year(ia)', mean(Q(ia,:), 2), '-', 'Color',[0.5 0.5 0.5], 'LineWidth',1.5)
plot(T.Year(ib), T.HEMO(ib), 'k-', 'LineWidth',1.5)
set(gca, 'TickDir','out', 'XTickLabel','')
box off;
ylabel('Q (m^{3} s^{-1})')
ylim = get(gca, 'YLim');
text(1931, ylim(2), 'A', 'FontSize',12, 'VerticalAlignment','middle', 'FontWeight','bold')
ax = gca; ax.Position(2) = ax.Position(2)+0.02;

subplot(2,4,4)
lim = [0 10500];
plot(lim, lim, '-', 'Color',[0.5 0.5 0.5])
hold on;
scatter(T.HEMO(ib), mean(Q(ia,:), 2), 10, "black","filled")
set(gca, 'TickDir','out', 'YLim',lim, 'XLim',lim, 'TickLength',[0.04 0])
box off;
ax = gca;
ax.Position(1) = 0.82;
ax.Position(2) = ax.Position(2)+0.02;
ax.Position(3) = 0.12;
ylabel('MsTMIP (m^{3} s^{-1})')
r = corr(T.HEMO(ib), mean(Q(ia,:), 2), 'rows','pairwise');
text(lim(2), lim(1), sprintf('R = %.02f', r), 'HorizontalAlignment','right', 'VerticalAlignment','bottom', 'FontSize',7)
text(diff(lim)*0.04, lim(2), 'B', 'FontSize',12, 'VerticalAlignment','middle', 'FontWeight','bold')

load ./data/McCabeWilliams_WaterBudget.mat;
[~, ia, ib] = intersect(year, T.Year);
Q = Q_S3_WY;

subplot(2,4,[5 7])
p1 = plot(year(ia)', Q(ia), '-', 'Color',[0.5 0.5 0.5], 'LineWidth',1.5);
hold on;
plot(T.Year(ib), T.HEMO(ib), 'k-', 'LineWidth',1.5)
set(gca, 'TickDir','out')
box off;
ylabel('Q (m^{3} s^{-1})')
xlabel('Water year')
ylim = get(gca, 'YLim');
text(1931, ylim(2), 'C', 'FontSize',12, 'VerticalAlignment','middle', 'FontWeight','bold')
ax = gca; ax.Position(2) = ax.Position(2)+0.02;

subplot(2,4,8)
lim = [0 8000];
plot(lim, lim, '-', 'Color',[0.5 0.5 0.5])
hold on;
scatter(T.HEMO(ib), mean(Q(ia,:), 2), 10, "black","filled")
set(gca, 'TickDir','out', 'YLim',lim, 'XLim',lim, 'TickLength',[0.04 0])
box off;
ax = gca;
ax.Position(1) = 0.82;
ax.Position(2) = ax.Position(2)+0.02;
ax.Position(3) = 0.12;
ylabel('MW11 (m^{3} s^{-1})')
xl = xlabel('Gage (m^{3} s^{-1})');
xl.Position(2) = -1500;
r = corr(T.HEMO(ib), mean(Q(ia,:), 2), 'rows','pairwise');
text(lim(2), lim(1), sprintf('R = %.02f', r), 'HorizontalAlignment','right', 'VerticalAlignment','bottom', 'FontSize',7)
text(diff(lim)*0.04, lim(2), 'D', 'FontSize',12, 'VerticalAlignment','middle', 'FontWeight','bold')

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r600','./output/models-vs-gages.tif')
close all;


