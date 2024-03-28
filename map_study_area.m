latlim=[36 50];
lonlim=[-115 -89];

states = shaperead('usastatehi','UseGeoCoords', true);
worldland= shaperead('landareas', 'UseGeoCoords', true);
load('./data/Mizu_streams_Strahler3_7.mat')
load('./data/MRB_subregions_GP.mat')

[x1, y1] = centroid(polyshape(SR1(1).Lon, SR1(1).Lat));
[x2, y2] = centroid(polyshape(SR2(2).Lon, SR2(2).Lat));
[x3, y3] = centroid(polyshape(SR3(2).Lon, SR3(2).Lat));
[x4, y4] = centroid(polyshape(SR4(2).Lon, SR4(2).Lat));
[x5, y5] = centroid(polyshape(SR5(2).Lon, SR5(2).Lat));
[x6, y6] = centroid(polyshape(SR6(2).Lon, SR6(2).Lat));

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 4.5 5];
clr = wesanderson('fantasticfox1');

% hydrology
subplot(3,1,[1 2])
axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',4,'MLineLocation',6,'MeridianLabel','on',...
        'ParallelLabel','on','GLineWidth',0.3,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','MLabelParallel',49.9999999,...
        'FLineWidth',1, 'GColor',[0.5 0.5 0.5], 'FontSize',7, 'FontColor',[0.5 0.5 0.5])
axis off;
axis image;
geoshow(SR1(1,1),'FaceColor',[0.8 0.8 0.8],'EdgeColor','none','LineWidth',1.2)
geoshow(SR2(2,1),'FaceColor',[0.8 0.8 0.8],'EdgeColor','none','LineWidth',1.2)
geoshow(SR3(2,1),'FaceColor',[0.8 0.8 0.8],'EdgeColor','none','LineWidth',1.2)
geoshow(SR4(2,1),'FaceColor',[0.8 0.8 0.8],'EdgeColor','none','LineWidth',1.2)
geoshow(SR5(2,1),'FaceColor',[0.8 0.8 0.8],'EdgeColor','none','LineWidth',1.2)
geoshow(SR6(2,1),'FaceColor',[0.8 0.8 0.8],'EdgeColor','none','LineWidth',1.2)
geoshow(states,'FaceColor','none','EdgeColor',[0.5 0.5 0.5],'LineWidth',0.5)

geoshow(MizuStrahler4,'Color',[33,102,172]/255,'LineWidth',0.5)
geoshow(MizuStrahler5,'Color',[33,102,172]/255,'LineWidth',1)
geoshow(MizuStrahler6,'Color',[33,102,172]/255,'LineWidth',1.25)
geoshow(MizuStrahler7,'Color',[33,102,172]/255,'LineWidth',1.5)
geoshow(SR1(1,1),'FaceColor','none','EdgeColor','k','LineWidth',1.)
geoshow(SR2(2,1),'FaceColor','none','EdgeColor','k','LineWidth',1.)
geoshow(SR3(2,1),'FaceColor','none','EdgeColor','k','LineWidth',1.)
geoshow(SR4(2,1),'FaceColor','none','EdgeColor','k','LineWidth',1.)
geoshow(SR5(2,1),'FaceColor','none','EdgeColor','k','LineWidth',1.)
geoshow(SR6(2,1),'FaceColor','none','EdgeColor','k','LineWidth',1.)
textm(y1, x1, '1', 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontWeight','bold', 'FontSize',16)
textm(y2, x2, '2', 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontWeight','bold', 'FontSize',16)
textm(y3, x3, '3', 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontWeight','bold', 'FontSize',16)
textm(y4, x4, '4', 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontWeight','bold', 'FontSize',16)
textm(y5, x5, '5', 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontWeight','bold', 'FontSize',16)
textm(y6, x6, '6', 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontWeight','bold', 'FontSize',16)

lat = [43.0529 42.4869 40.6818 39.1123 38.9789 38.7094];
lon = [-98.5621 -96.4122 -95.846 -94.5877 -92.754 -91.4379];
gage = {'Fort Randall Dam', 'Sioux City, IA', 'Nebraska City, NE', 'Kansas City, MO', 'Boonville, MO', 'Hermann, MO'};

scatterm(lat(6), lon(6), 50, "black", "filled", 'MarkerFaceColor','w', 'MarkerEdgeColor','k', 'Marker','o', 'LineWidth',2)
textm(lat(6)+0.5, lon(6)-0.9, 'Hermann, MO','HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',8)

ax = gca;
subplotsqueeze(ax, 1.15)
text(-0.169,0.84,'a', 'FontSize',12)
ax.Position(2) = 0.4;

% gage
subplot(3,1,3)
T = readtable('./data/Missouri_Big6_UpdatedUnregulatedFlow_DailyMonthlyWY_avgCMS.xlsx',...
    'Sheet','WaterYear_cms');
T.Properties.VariableNames = {'Year','FTRA','SUX','NCNE','MKC','MNMO','HEMO'};
mdl = fitlm([1931:2010]', T.HEMO(T.Year>=1931 & T.Year<=2010));
beta(1) = mdl.Coefficients.Estimate(2);
se(1) = mdl.Coefficients.SE(2);

% MsTMIP
load ./data/MsTMIP_WaterBudget;
[~, ia, ib] = intersect(year, T.Year);
Q = WB_SG3_unscaled;
dQ = Q - repmat(mean(Q(year>=1931&year<=1960, :)), length(year), 1);
ci = quantile(dQ, [0.05 0.95],2);
r_mstmip = corr(mean(dQ(ia,:),2), T.HEMO(ib));
mdl = fitlm([1931:2010]', mean(dQ(year>=1931 & year<=2010,:), 2));
beta(2) = mdl.Coefficients.Estimate(2);
se(2) = mdl.Coefficients.SE(2);

fill([year(2:end) fliplr(year(2:end))],...
    [ci(2:end,1)' fliplr(ci(2:end,2)')],...
    clr(1,:), 'EdgeColor','none', 'FaceAlpha',0.3)
hold on;
p2 = plot(year', mean(dQ, 2), '-', 'Color',clr(1,:), 'LineWidth',1.5);

% WBM
load ./data/McCabeWilliams_WaterBudget.mat;
[~, ia, ib] = intersect(year, T.Year);
Q = Q_S3_WY;
dQ = Q-mean(Q(year>=1931 & year<=1960));
r_wbm = corr(dQ(ia), T.HEMO(ib));
mdl = fitlm([1931:2010]', dQ(year>=1931 & year<=2010));
beta(3) = mdl.Coefficients.Estimate(2);
se(3) = mdl.Coefficients.SE(2);
p3 = plot(year', dQ, '-', 'Color',clr(2,:), 'LineWidth',1.5);
p1 = plot(T.Year, T.HEMO-mean(T.HEMO(T.Year>=1931 & T.Year<=1960)), 'k-', 'LineWidth',1.5);

set(gca, 'TickDir','out', 'XLim', [1930 2020], 'TickLength',[0.015 0.],'YLim',[-4000 8000])
box off;
yl = ylabel('\DeltaQ (m^{3} s^{-1})', 'FontSize',11);
yl.Position(1) = 1919;
ax = gca;
ax.Position(1) = 0.17;
ax.Position(2) = 0.08;
ax.Position(3) = 0.775;
ax.Position(4) = 0.23;
text(1916,ax.YLim(2),'b', 'FontSize',12, 'VerticalAlignment','baseline')

lgd = legend([p1 p2 p3],...
    sprintf('gage (trend = %0.1f \x00B1 %0.1f m^3 s^{-1} yr^{-1})',beta(1),se(1)*1.96),...
    sprintf('MsTMIP (R = %0.02f; trend = %0.1f \x00B1 %0.1f m^3 s^{-1} yr^{-1})',r_mstmip,beta(2),se(2)*1.96),...
    sprintf('MW11 (R = %0.02f; trend = %0.1f \x00B1 %0.1f m^3 s^{-1} yr^{-1})',r_wbm,beta(3),se(3)*1.96), 'Location','northwest');
legend('boxoff')
lgd.Position(2) = 0.28;

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/mrb-study-area.tif')
close all;


