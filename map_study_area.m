latlim=[36 50];
lonlim=[-115 -89];

states = shaperead('usastatehi','UseGeoCoords', true);
worldland= shaperead('landareas', 'UseGeoCoords', true);
load('./data/Mizu_streams_Strahler3_7.mat')
load('./data/MRB_subregions_GP.mat')

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 4.5 5];

% hydrology
subplot(3,1,[1 2])
axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',4,'MLineLocation',6,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.3,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','MLabelParallel','north',...
        'FLineWidth',1, 'GColor',[0.5 0.5 0.5], 'FontSize',8)
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

lat = [43.0529 42.4869 40.6818 39.1123 38.9789 38.7094];
lon = [-98.5621 -96.4122 -95.846 -94.5877 -92.754 -91.4379];
gage = {'Fort Randall Dam', 'Sioux City, IA', 'Nebraska City, NE', 'Kansas City, MO', 'Boonville, MO', 'Hermann, MO'};

scatterm(lat(1:5), lon(1:5), 40, "black", "filled", 'MarkerFaceColor','w', 'MarkerEdgeColor','k', 'Marker','o')
scatterm(lat(6), lon(6), 150, "black", "filled", 'MarkerFaceColor','w', 'MarkerEdgeColor','k', 'Marker','pentagram')

ax = gca;
subplotsqueeze(ax, 1.2)
text(-0.16,0.83,'a', 'FontSize',12)

subplot(3,1,3)
T = readtable('./data/Missouri_Big6_UpdatedUnregulatedFlow_DailyMonthlyWY_avgCMS.xlsx',...
    'Sheet','WaterYear_cms');
T.Properties.VariableNames = {'Year','FTRA','SUX','NCNE','MKC','MNMO','HEMO'};
plot(T.Year, T.HEMO, 'k-', 'LineWidth',1.5)
hold on;
plot([1931 2019], [mean(T.HEMO(T.Year>=1931 & T.Year<=1960)) mean(T.HEMO(T.Year>=1931 & T.Year<=1960))], 'k-', 'LineWidth',0.5)
plot([1981 2019], [mean(T.HEMO(T.Year>=1981 & T.Year<=2019)) mean(T.HEMO(T.Year>=1981 & T.Year<=2019))], 'k--', 'LineWidth',0.75)
hold off;
set(gca, 'TickDir','out', 'XLim', [1930 2020], 'TickLength',[0.015 0.])
box off;
yl = ylabel('Q (m^{3} s^{-1})', 'FontSize',11);
yl.Position(1) = 1919;
ax = gca;
ax.Position(1) = 0.17;
ax.Position(2) = 0.08;
ax.Position(3) = 0.775;
ax.Position(4) = 0.23;
text(1916.5,ax.YLim(2),'b', 'FontSize',12)

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/mrb-study-area.tif')
close all;


