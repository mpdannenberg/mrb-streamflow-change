latlim=[36 50];
lonlim=[-115 -89];

states = shaperead('usastatehi','UseGeoCoords', true);
worldland= shaperead('landareas', 'UseGeoCoords', true);
load('./data/Mizu_streams_Strahler3_7.mat')
load('./data/MRB_subregions_GP.mat')

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 3.5 5.5];

% hydrology
subplot(5,1,[1 2])
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

lat = [43.0529 42.4869 40.6818 39.1123 38.9789 38.7094];
lon = [-98.5621 -96.4122 -95.846 -94.5877 -92.754 -91.4379];
gage = {'Fort Randall Dam', 'Sioux City, IA', 'Nebraska City, NE', 'Kansas City, MO', 'Boonville, MO', 'Hermann, MO'};

scatterm(lat(1:5), lon(1:5), 20, "black", "filled", 'MarkerFaceColor','w', 'MarkerEdgeColor','k', 'Marker','o')
scatterm(lat(6), lon(6), 75, "black", "filled", 'MarkerFaceColor','w', 'MarkerEdgeColor','k', 'Marker','pentagram')

ax = gca;
ax.Position(2) = 0.68;
text(-0.16,0.83,'a', 'FontSize',12)

% land cover
[dat, R] = readgeoraster("data\nlcd_sub2.tif");
T = table('Size',[15 2], 'VariableTypes',{'double','double'}, 'VariableNames',{'Class','Proportion'});
T.Class = [11 21 22 23 24 31 41 42 43 52 71 81 82 90 95]';
N = sum(sum(dat>0 & dat<100));
writetable(T, './output/nlcd-proportion.xlsx');

dat(dat>100) = 0;
[ny,nx] = size(dat);
temp = double(dat(1:4:ny, 1:4:nx));
temp(temp==0) = NaN;
for i = 1:size(T, 1)
    temp(temp == T.Class(i)) = i;
end
cellsize = R.CellExtentInLatitude*4;
lat = (R.LatitudeLimits(2)-cellsize/2):(-cellsize):(R.LatitudeLimits(1)+cellsize/2);
lon = (R.LongitudeLimits(1)+cellsize/2):(cellsize):(R.LongitudeLimits(2)-cellsize/2);

%dat(dat > 100) = NaN;
clr = [71 107 161
    222 201 201
    217 148 130
    237 0 0
    171 0 0
    179 173 163
    105 171 99
    28 99 48
    181 201 143
    204 186 125
    227 227 194
    219 217 61
    171 112 41
    186 217 235
    112 163 186]/255;
subplot(5,1,[3 4])
axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',4,'MLineLocation',6,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.3,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','MLabelParallel','north',...
        'FLineWidth',1, 'GColor',[0.5 0.5 0.5], 'FontSize',8)
axis off;
axis image;
surfm(lat, lon, temp)
geoshow(states,'FaceColor','none','EdgeColor',[0.5 0.5 0.5],'LineWidth',0.5)
geoshow(SR1(1,1),'FaceColor','none','EdgeColor','k','LineWidth',1.2)
geoshow(SR2(2,1),'FaceColor','none','EdgeColor','k','LineWidth',1.2)
geoshow(SR3(2,1),'FaceColor','none','EdgeColor','k','LineWidth',1.2)
geoshow(SR4(2,1),'FaceColor','none','EdgeColor','k','LineWidth',1.2)
geoshow(SR5(2,1),'FaceColor','none','EdgeColor','k','LineWidth',1.2)
geoshow(SR6(2,1),'FaceColor','none','EdgeColor','k','LineWidth',1.2)
caxis([0.5 15.5])
colormap(clr)
ax = gca;
ax.Position(2) = 0.35;
text(-0.16,0.83,'b', 'FontSize',12)

cb = colorbar('southoutside');
cb.Position(1) = 0.17;
cb.Position(2) = 0.3;
cb.Position(3) = 0.7;
cb.FontSize = 8;
cb.Ticks = 1:15;
cb.TickLength = 0;
cb.TickLabels = {'11','21','22','23','24','31','41','42','43','52','71','81','82','90','95'};
xlabel(cb, 'NLCD cover class')

subplot(5,1,5)
T = readtable('./data/Missouri_Big6_UpdatedUnregulatedFlow_DailyMonthlyWY_avgCMS.xlsx',...
    'Sheet','WaterYear_cms');
T.Properties.VariableNames = {'Year','FTRA','SUX','NCNE','MKC','MNMO','HEMO'};
plot(T.Year, T.HEMO, 'k-', 'LineWidth',1.5)
hold on;
plot([1931 2019], [mean(T.HEMO(T.Year>=1931 & T.Year<=1960)) mean(T.HEMO(T.Year>=1931 & T.Year<=1960))], 'k-', 'LineWidth',0.5)
plot([1981 2019], [mean(T.HEMO(T.Year>=1981 & T.Year<=2019)) mean(T.HEMO(T.Year>=1981 & T.Year<=2019))], 'k--', 'LineWidth',0.75)
hold off;
set(gca, 'TickDir','out', 'XLim', [1930 2020])
box off;
ylabel('Q (m^{3} s^{-1})')
ax = gca;
ax.Position(1) = 0.16;
ax.Position(2) = 0.08;
ax.Position(3) = 0.72;
text(1934,6000,'c', 'FontSize',12)

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/mrb-study-area.tif')
close all;


