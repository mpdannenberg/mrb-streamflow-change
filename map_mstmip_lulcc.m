% Map change in MsTMIP land cover

%% Basic script parameters
syear = 1801;
eyear = 2010;
load ./data/MRB_subregions_GP;
latlim = [32 50];
lonlim = [-115 -85];

%% Get MsTMIP LULCC data
cd('D:/Data_Analysis/MsTMIP/')
info = ncinfo('mstmip_driver_global_hd_lulcc_1801_v1.nc4');
lat = double(ncread('mstmip_driver_global_hd_lulcc_1801_v1.nc4', 'lat')); latidx = lat >= latlim(1) & lat <= latlim(2); ny = sum(latidx);
lon = double(ncread('mstmip_driver_global_hd_lulcc_1801_v1.nc4', 'lon')); lonidx = lon >= lonlim(1) & lon <= lonlim(2); nx = sum(lonidx);

%% Find MRB grid cells and calculate area of each cell
[LON, LAT] = meshgrid(lon(lonidx), lat(latidx));
LatLon = [reshape(LAT, [], 1) reshape(LON, [], 1)];
idx = zeros(size(LatLon,1),1);
[IN, ON] = inpolygon(LatLon(:,1), LatLon(:,2), SR6(2).Lat, SR6(2).Lon);
idx = idx | IN | ON;
[IN, ON] = inpolygon(LatLon(:,1), LatLon(:,2), SR1(1).Lat, SR1(1).Lon);
idx = idx | IN | ON;
[IN, ON] = inpolygon(LatLon(:,1), LatLon(:,2), SR2(2).Lat, SR2(2).Lon);
idx = idx | IN | ON;
[IN, ON] = inpolygon(LatLon(:,1), LatLon(:,2), SR3(2).Lat, SR3(2).Lon);
idx = idx | IN | ON;
[IN, ON] = inpolygon(LatLon(:,1), LatLon(:,2), SR4(2).Lat, SR4(2).Lon);
idx = idx | IN | ON;
[IN, ON] = inpolygon(LatLon(:,1), LatLon(:,2), SR5(2).Lat, SR5(2).Lon);
idx = idx | IN | ON;
MRBidx = reshape(idx, ny, nx);

e = referenceEllipsoid('World Geodetic System 1984');
area = areaquad(reshape(LAT-0.25,[],1),reshape(LON-0.25,[],1),reshape(LAT+0.25,[],1),reshape(LON+0.25,[],1),e);
area = reshape(area, ny, nx); 

clear LAT LON LatLon IN ON idx e;

%% Get LULCC data at beginning and end
type = ncread('mstmip_driver_global_hd_lulcc_1801_v1.nc4', 'type');
lc1801 = permute(ncread('mstmip_driver_global_hd_lulcc_1801_v1.nc4', 'biome_frac'), [2 1 3]);
tree1801 = sum(lc1801(latidx, lonidx, type>=1 & type<=9), 3) + 0.5*sum(lc1801(latidx, lonidx, type>=10 & type<=36), 3);
shrub1801 = lc1801(latidx, lonidx, type==37) + 0.5*sum(lc1801(latidx, lonidx, (type>=10 & type<=18) | (type>=38 & type<=40)), 3);
grass1801 = lc1801(latidx, lonidx, type==41) + 0.5*sum(lc1801(latidx, lonidx, (type>=19 & type<=27) | type==38 | (type>=42 & type<=43)), 3);
crop1801 = lc1801(latidx, lonidx, type==44) + 0.5*sum(lc1801(latidx, lonidx, (type>=28 & type<=36) | type==39 | type==43), 3);
tree1801(~MRBidx) = NaN; shrub1801(~MRBidx) = NaN; grass1801(~MRBidx) = NaN; crop1801(~MRBidx) = NaN; 

lc2010 = permute(ncread('mstmip_driver_global_hd_lulcc_2010_v1.nc4', 'biome_frac'), [2 1 3]);
tree2010 = sum(lc2010(latidx, lonidx, type>=1 & type<=9), 3) + 0.5*sum(lc2010(latidx, lonidx, type>=10 & type<=36), 3);
shrub2010 = lc2010(latidx, lonidx, type==37) + 0.5*sum(lc2010(latidx, lonidx, (type>=10 & type<=18) | (type>=38 & type<=40)), 3);
grass2010 = lc2010(latidx, lonidx, type==41) + 0.5*sum(lc2010(latidx, lonidx, (type>=19 & type<=27) | type==38 | (type>=42 & type<=43)), 3);
crop2010 = lc2010(latidx, lonidx, type==44) + 0.5*sum(lc2010(latidx, lonidx, (type>=28 & type<=36) | type==39 | type==43), 3);
tree2010(~MRBidx) = NaN; shrub2010(~MRBidx) = NaN; grass2010(~MRBidx) = NaN; crop2010(~MRBidx) = NaN; 

%% Get LULCC time series
yr = syear:eyear;
nt = length(yr);
tree = NaN(size(yr));
shrub = NaN(size(yr));
crop = NaN(size(yr));
grass = NaN(size(yr));

totalarea = sum(area(MRBidx));

for i = 1:nt
    
    lc = permute(ncread(['mstmip_driver_global_hd_lulcc_',num2str(yr(i)),'_v1.nc4'], 'biome_frac'), [2 1 3]);
    t = sum(lc(latidx, lonidx, type>=1 & type<=9), 3) + 0.5*sum(lc(latidx, lonidx, type>=10 & type<=36), 3);
    s = lc(latidx, lonidx, type==37) + 0.5*sum(lc(latidx, lonidx, (type>=10 & type<=18) | (type>=38 & type<=40)), 3);
    g = lc(latidx, lonidx, type==41) + 0.5*sum(lc(latidx, lonidx, (type>=19 & type<=27) | type==38 | (type>=42 & type<=43)), 3);
    c = lc(latidx, lonidx, type==44) + 0.5*sum(lc(latidx, lonidx, (type>=28 & type<=36) | type==39 | type==43), 3);

    tree(i) = sum(t(MRBidx) .* area(MRBidx)) / totalarea;
    shrub(i) = sum(s(MRBidx) .* area(MRBidx)) / totalarea;
    crop(i) = sum(c(MRBidx) .* area(MRBidx)) / totalarea;
    grass(i) = sum(g(MRBidx) .* area(MRBidx)) / totalarea;
    
end


%% Map a map of LULCC change (1801-2010)
cd('D:\Publications\Dannenberg_et_al_MRB_streamflow_change')
states = shaperead('usastatehi','UseGeoCoords', true);
latlim=[36 50];
lonlim=[-115 -89];

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6.5 6.5];
clr = cbrewer('div','RdBu',20); clr(clr<0) = 0; clr(clr>1) = 1;

subplot(3,2,1)
axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',4,'MLineLocation',6,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.3,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','MLabelParallel','north',...
        'FLineWidth',1, 'GColor',[0.5 0.5 0.5], 'FontSize',8)
axis off;
axis image;
surfm(double(lat(latidx))+0.25, double(lon(lonidx))-0.25, tree2010-tree1801);
caxis([-1 1])
colormap(gca, clr)
geoshow(states,'FaceColor','none','EdgeColor',[0.5 0.5 0.5],'LineWidth',0.5)
geoshow(SR1(1,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR2(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR3(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR4(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR5(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR6(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
ttl = title('tree', 'FontSize',12); ttl.Position(2) = 0.84;
ax = gca;
ax.Position(1) = 0.1;
ax.Position(2) = 0.7;
subplotsqueeze(ax, 1.1)
text(-0.16, 0.83, 'A', 'FontSize',12, 'FontWeight','bold')

subplot(3,2,2)
axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',4,'MLineLocation',6,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.3,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','MLabelParallel','north',...
        'FLineWidth',1, 'GColor',[0.5 0.5 0.5], 'FontSize',8)
axis off;
axis image;
surfm(double(lat(latidx))+0.25, double(lon(lonidx))-0.25, shrub2010-shrub1801);
caxis([-1 1])
colormap(gca, clr)
geoshow(states,'FaceColor','none','EdgeColor',[0.5 0.5 0.5],'LineWidth',0.5)
geoshow(SR1(1,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR2(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR3(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR4(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR5(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR6(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
ttl = title('shrub', 'FontSize',12); ttl.Position(2) = 0.84;
ax = gca;
ax.Position(1) = 0.56;
ax.Position(2) = 0.7;
subplotsqueeze(ax, 1.1)
text(-0.16, 0.83, 'B', 'FontSize',12, 'FontWeight','bold')

subplot(3,2,3)
axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',4,'MLineLocation',6,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.3,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','MLabelParallel','north',...
        'FLineWidth',1, 'GColor',[0.5 0.5 0.5], 'FontSize',8)
axis off;
axis image;
surfm(double(lat(latidx))+0.25, double(lon(lonidx))-0.25, grass2010-grass1801);
caxis([-1 1])
colormap(gca, clr)
geoshow(states,'FaceColor','none','EdgeColor',[0.5 0.5 0.5],'LineWidth',0.5)
geoshow(SR1(1,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR2(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR3(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR4(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR5(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR6(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
ttl = title('grass', 'FontSize',12); ttl.Position(2) = 0.84;
ax = gca;
ax.Position(1) = 0.1;
ax.Position(2) = 0.4;
subplotsqueeze(ax, 1.1)
text(-0.16, 0.83, 'C', 'FontSize',12, 'FontWeight','bold')

subplot(3,2,4)
axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',4,'MLineLocation',6,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.3,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','MLabelParallel','north',...
        'FLineWidth',1, 'GColor',[0.5 0.5 0.5], 'FontSize',8)
axis off;
axis image;
surfm(double(lat(latidx))+0.25, double(lon(lonidx))-0.25, crop2010-crop1801);
caxis([-1 1])
colormap(gca, clr)
geoshow(states,'FaceColor','none','EdgeColor',[0.5 0.5 0.5],'LineWidth',0.5)
geoshow(SR1(1,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR2(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR3(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR4(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR5(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR6(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
ttl = title('crop', 'FontSize',12); ttl.Position(2) = 0.84;
ax = gca;
ax.Position(1) = 0.56;
ax.Position(2) = 0.4;
subplotsqueeze(ax, 1.1)
text(-0.16, 0.83, 'D', 'FontSize',12, 'FontWeight','bold')

cb = colorbar('southoutside');
cb.Position = [0.1 0.34 0.8 0.025];
cb.Ticks = -1:0.1:1;
cb.TickLabels = {'-1','','-0.8','','-0.6','','-0.4','','-0.2','','0','','0.2','','0.4','','0.6','','0.8','','1'};
cb.TickLength = 0.03;
xlabel(cb, '\Deltaf_{c}', 'FontSize',10)

%% time series
subplot(3,2,[5 6])
ax = gca;
ax.Position(2) = 0.06;
ax.Position(4) = 0.18;

clr = [105 171 99
    28 99 48
    204 186 125
    171 112 41]/255;
plot(yr, tree, '-','Color',clr(2,:), 'LineWidth',2)
hold on;
plot(yr, grass, '-','Color',clr(1,:), 'LineWidth',2)
plot(yr, crop, '-','Color',clr(4,:), 'LineWidth',2)
plot(yr, shrub, '-','Color',clr(3,:), 'LineWidth',2)
set(ax, 'XLim',[syear eyear], 'TickDir','out')
box off;
lgd = legend('tree','grass','crop','shrub', 'Location','north', 'Orientation','horizontal');
lgd.Position(2) = 0.22;
legend('boxoff')
yl = ylabel('f_{c}', 'Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle', 'FontSize',10);
yl.Position(1) = 1785;
xlabel('Year')
text(1805, 0.6, 'E', 'FontSize',12, 'FontWeight','bold')
hold off;

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/MsTMIP-lulcc-mrb.tif')
close all;

