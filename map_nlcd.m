% map NLCD land cover

latlim=[36 50];
lonlim=[-115 -89];

states = shaperead('usastatehi','UseGeoCoords', true);
worldland= shaperead('landareas', 'UseGeoCoords', true);
load('./data/MRB_subregions_GP.mat')

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6.5 5];

% calculate areal proportion
[dat, R] = readgeoraster("data\nlcd_sub2.tif");
T = table('Size',[15 2], 'VariableTypes',{'double','double'}, 'VariableNames',{'Class','Proportion'});
T.Class = [11 21 22 23 24 31 41 42 43 52 71 81 82 90 95]';
N = sum(sum(dat>0 & dat<100));
for i = 1:size(T, 1)
    T.Proportion(i) = sum(sum(dat == T.Class(i))) / N;
end
writetable(T, './output/nlcd-proportion.xlsx');

% reclass and sample every fourth to reduce size for mapping
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

% map
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
axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',4,'MLineLocation',6,'MeridianLabel','on',...
        'ParallelLabel','on','GLineWidth',0.3,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','MLabelParallel',49.99999999,...
        'FLineWidth',1, 'GColor',[0.5 0.5 0.5], 'FontSize',7, 'FontColor',[0.5 0.5 0.5])
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
ax.Position(2) = 0.22;

cb = colorbar('southoutside');
cb.Position(1) = 0.13;
cb.Position(2) = 0.2;
cb.Position(3) = 0.78;
cb.FontSize = 8;
cb.Ticks = 1:15;
cb.TickLength = 0;
cb.TickLabels = {'Open water','Developed, open space','Developed, low intensity',...
    'Developed, medium intensity','Developed, high intensity','Barren land',...
    'Deciduous forest','Evergreen forest','Mixed forest','Shrub/scrub',...
    'Grassland/herbaceous','Pasture/hay','Cultivated crops',...
    'Woody wetlands','Emergent herbaceous westlands'};
xs = cb.Position(1);
xe = cb.Position(3) + cb.Position(1);
xw = (xe - xs) / 15;
ys = cb.Position(2);
yw = cb.Position(4);
for i = 1:15
    x = xs + (xw/2) + xw*(i-1);
    annotation('textbox',[xs+xw*(i-1) ys xw yw],'String',sprintf('%0.1f',100*T.Proportion(i)),...
        'BackgroundColor','none', "FontSize",8, "VerticalAlignment","middle", "HorizontalAlignment","center",'LineWidth',1.1)
end
annotation('textbox',[0 ys xs-0.01 yw],'String','f_{c} (%)',...
        'BackgroundColor','none', "FontSize",10, "VerticalAlignment","middle",...
        "HorizontalAlignment","right",'Color','k', 'EdgeColor','none')

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/mrb-nlcd.tif')
close all;


