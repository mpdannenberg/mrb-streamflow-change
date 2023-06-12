
T = readtable('./data/Missouri_Big6_UpdatedUnregulatedFlow_DailyMonthlyWY_avgCMS.xlsx',...
    'Sheet','WaterYear_cms');
T.Properties.VariableNames = {'Year','FTRA','SUX','NCNE','MKC','MNMO','HEMO'};

basins = shaperead('./data/LMRBgagebasins/lmrb_shapefiles.shp');

%% Find WBM grid cells for each basin and calculate area
load ./data/McCabeWilliams_WaterBudget.mat;
[nt,ny,nx] = size(Rs_S3_WY);

[LON, LAT] = meshgrid(lon, lat);
LatLon = [reshape(LAT, [], 1) reshape(LON, [], 1)];

[IN, ON] = inpolygon(LatLon(:,1), LatLon(:,2), basins(1).Y, basins(1).X);
idx = IN | ON;
FTRA = reshape(idx, ny, nx);

[IN, ON] = inpolygon(LatLon(:,1), LatLon(:,2), basins(2).Y, basins(2).X);
idx = IN | ON;
SUX = reshape(idx, ny, nx);

[IN, ON] = inpolygon(LatLon(:,1), LatLon(:,2), basins(3).Y, basins(3).X);
idx = IN | ON;
NCNE = reshape(idx, ny, nx);

[IN, ON] = inpolygon(LatLon(:,1), LatLon(:,2), basins(4).Y, basins(4).X);
idx = IN | ON;
MKC = reshape(idx, ny, nx);

[IN, ON] = inpolygon(LatLon(:,1), LatLon(:,2), basins(5).Y, basins(5).X);
idx = IN | ON;
MNMO = reshape(idx, ny, nx);

[IN, ON] = inpolygon(LatLon(:,1), LatLon(:,2), basins(6).Y, basins(6).X);
idx = IN | ON;
HEMO = reshape(idx, ny, nx);

e = referenceEllipsoid('World Geodetic System 1984');
area = areaquad(reshape(LAT-0.125,[],1),reshape(LON-0.125,[],1),reshape(LAT+0.125,[],1),reshape(LON+0.125,[],1),e);
area = reshape(area, ny, nx); 

clear LAT LON LatLon IN ON idx e;

%% Make figure showing comparison between gage data and WBM
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6.5 6];

[~, ia, ib] = intersect(year, T.Year);

% FTRA
Q = sum((0.001 .* Rs_S3_WY(:, FTRA) .* repmat(area(FTRA)', size(Rs_S3_WY, 1), 1)) / (365 * 24 * 60 * 60), 2);

subplot(6,4,[1 3])
p1 = plot(year(ia)', Q(ia), '-', 'Color',[0.5 0.5 0.5], 'LineWidth',1.5);
hold on;
p2 = plot(T.Year(ib), T.FTRA(ib), 'k-', 'LineWidth',1.5);
set(gca, 'TickDir','out', 'XTickLabel','')
box off;
ylabel('Q (m^{3} s^{-1})')
ylim = get(gca, 'YLim');
text(1931, ylim(2), 'Fort Randall Dam', 'FontSize',10, 'VerticalAlignment','middle')
lgd = legend([p2 p1], 'Gage', 'WBM', 'Location','north', 'Orientation','horizontal');
lgd.Position(1) = 0.45;
lgd.Position(2) = 0.92;
legend('boxoff')

subplot(6,4,4)
lim = [0 2800];
plot(lim, lim, '-', 'Color',[0.5 0.5 0.5])
hold on;
scatter(T.FTRA(ib), mean(Q(ia,:), 2), 10, "black","filled")
set(gca, 'TickDir','out', 'YLim',lim, 'XLim',lim, 'TickLength',[0.04 0])
box off;
ax = gca;
ax.Position(1) = 0.82;
ax.Position(3) = 0.12;
ylabel('WBM')
r = corr(T.FTRA(ib), mean(Q(ia,:), 2), 'rows','pairwise');
text(lim(2), lim(1), sprintf('R = %.02f', r), 'HorizontalAlignment','right', 'VerticalAlignment','bottom', 'FontSize',7)

% SUX
Q = sum((0.001 .* Rs_S3_WY(:, SUX) .* repmat(area(SUX)', size(Rs_S3_WY, 1), 1)) / (365 * 24 * 60 * 60), 2);

subplot(6,4,[5 7])
p1 = plot(year(ia)', Q(ia), '-', 'Color',[0.5 0.5 0.5], 'LineWidth',1.5);
hold on;
plot(T.Year(ib), T.SUX(ib), 'k-', 'LineWidth',1.5)
set(gca, 'TickDir','out', 'XTickLabel','')
box off;
ylabel('Q (m^{3} s^{-1})')
ylim = get(gca, 'YLim');
text(1931, ylim(2), 'Sioux City, IA', 'FontSize',10, 'VerticalAlignment','middle')

subplot(6,4,8)
lim = [0 3000];
plot(lim, lim, '-', 'Color',[0.5 0.5 0.5])
hold on;
scatter(T.SUX(ib), mean(Q(ia,:), 2), 10, "black","filled")
set(gca, 'TickDir','out', 'YLim',lim, 'XLim',lim, 'TickLength',[0.04 0])
box off;
ax = gca;
ax.Position(1) = 0.82;
ax.Position(3) = 0.12;
ylabel('WBM')
r = corr(T.SUX(ib), mean(Q(ia,:), 2), 'rows','pairwise');
text(lim(2), lim(1), sprintf('R = %.02f', r), 'HorizontalAlignment','right', 'VerticalAlignment','bottom', 'FontSize',7)

% NCNE
Q = sum((0.001 .* Rs_S3_WY(:, NCNE) .* repmat(area(NCNE)', size(Rs_S3_WY, 1), 1)) / (365 * 24 * 60 * 60), 2, "omitnan");

subplot(6,4,[9 11])
p1 = plot(year(ia)', Q(ia), '-', 'Color',[0.5 0.5 0.5], 'LineWidth',1.5);
hold on;
plot(T.Year(ib), T.NCNE(ib), 'k-', 'LineWidth',1.5)
set(gca, 'TickDir','out', 'XTickLabel','')
box off;
ylabel('Q (m^{3} s^{-1})')
ylim = get(gca, 'YLim');
text(1931, ylim(2), 'Nebraska City, NE', 'FontSize',10, 'VerticalAlignment','middle')

subplot(6,4,12)
lim = [0 4000];
plot(lim, lim, '-', 'Color',[0.5 0.5 0.5])
hold on;
scatter(T.NCNE(ib), mean(Q(ia,:), 2), 10, "black","filled")
set(gca, 'TickDir','out', 'YLim',lim, 'XLim',lim, 'TickLength',[0.04 0])
box off;
ax = gca;
ax.Position(1) = 0.82;
ax.Position(3) = 0.12;
ylabel('WBM')
r = corr(T.NCNE(ib), mean(Q(ia,:), 2), 'rows','pairwise');
text(lim(2), lim(1), sprintf('R = %.02f', r), 'HorizontalAlignment','right', 'VerticalAlignment','bottom', 'FontSize',7)

% MKC
Q = sum((0.001 .* Rs_S3_WY(:, MKC) .* repmat(area(MKC)', size(Rs_S3_WY, 1), 1)) / (365 * 24 * 60 * 60), 2, 'omitnan');

subplot(6,4,[13 15])
p1 = plot(year(ia)', Q(ia), '-', 'Color',[0.5 0.5 0.5], 'LineWidth',1.5);
hold on;
plot(T.Year(ib), T.MKC(ib), 'k-', 'LineWidth',1.5)
set(gca, 'TickDir','out', 'XTickLabel','')
box off;
ylabel('Q (m^{3} s^{-1})')
ylim = get(gca, 'YLim');
text(1931, ylim(2), 'Kansas City, MO', 'FontSize',10, 'VerticalAlignment','middle')

subplot(6,4,16)
lim = [0 6000];
plot(lim, lim, '-', 'Color',[0.5 0.5 0.5])
hold on;
scatter(T.MKC(ib), mean(Q(ia,:), 2), 10, "black","filled")
set(gca, 'TickDir','out', 'YLim',lim, 'XLim',lim, 'TickLength',[0.04 0])
box off;
ax = gca;
ax.Position(1) = 0.82;
ax.Position(3) = 0.12;
ylabel('WBM')
r = corr(T.MKC(ib), mean(Q(ia,:), 2), 'rows','pairwise');
text(lim(2), lim(1), sprintf('R = %.02f', r), 'HorizontalAlignment','right', 'VerticalAlignment','bottom', 'FontSize',7)

% MNMO
Q = sum((0.001 .* Rs_S3_WY(:, MNMO) .* repmat(area(MNMO)', size(Rs_S3_WY, 1), 1)) / (365 * 24 * 60 * 60), 2, 'omitnan');

subplot(6,4,[17 19])
p1 = plot(year(ia)', Q(ia), '-', 'Color',[0.5 0.5 0.5], 'LineWidth',1.5);
hold on;
plot(T.Year(ib), T.MNMO(ib), 'k-', 'LineWidth',1.5)
set(gca, 'TickDir','out', 'XTickLabel','')
box off;
ylabel('Q (m^{3} s^{-1})')
ylim = get(gca, 'YLim');
text(1931, ylim(2), 'Boonville, MO', 'FontSize',10, 'VerticalAlignment','middle')

subplot(6,4,20)
lim = [0 6000];
plot(lim, lim, '-', 'Color',[0.5 0.5 0.5])
hold on;
scatter(T.MNMO(ib), mean(Q(ia,:), 2), 10, "black","filled")
set(gca, 'TickDir','out', 'YLim',lim, 'XLim',lim, 'TickLength',[0.04 0])
box off;
ax = gca;
ax.Position(1) = 0.82;
ax.Position(3) = 0.12;
ylabel('WBM')
r = corr(T.MNMO(ib), mean(Q(ia,:), 2), 'rows','pairwise');
text(lim(2), lim(1), sprintf('R = %.02f', r), 'HorizontalAlignment','right', 'VerticalAlignment','bottom', 'FontSize',7)

% HEMO
Q = sum((0.001 .* Rs_S3_WY(:, HEMO) .* repmat(area(HEMO)', size(Rs_S3_WY, 1), 1)) / (365 * 24 * 60 * 60), 2, 'omitnan');

subplot(6,4,[21 23])
p1 = plot(year(ia)', Q(ia), '-', 'Color',[0.5 0.5 0.5], 'LineWidth',1.5);
hold on;
plot(T.Year(ib), T.HEMO(ib), 'k-', 'LineWidth',1.5)
set(gca, 'TickDir','out')
box off;
ylabel('Q (m^{3} s^{-1})')
ylim = get(gca, 'YLim');
text(1931, ylim(2), 'Hermann, MO', 'FontSize',10, 'VerticalAlignment','middle')

subplot(6,4,24)
lim = [0 8000];
plot(lim, lim, '-', 'Color',[0.5 0.5 0.5])
hold on;
scatter(T.HEMO(ib), mean(Q(ia,:), 2), 10, "black","filled")
set(gca, 'TickDir','out', 'YLim',lim, 'XLim',lim, 'TickLength',[0.04 0])
box off;
ax = gca;
ax.Position(1) = 0.82;
ax.Position(3) = 0.12;
ylabel('WBM')
xlabel('Gage')
r = corr(T.HEMO(ib), mean(Q(ia,:), 2), 'rows','pairwise');
text(lim(2), lim(1), sprintf('R = %.02f', r), 'HorizontalAlignment','right', 'VerticalAlignment','bottom', 'FontSize',7)

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r600','./output/WBM-vs-gages.tif')
close all;


