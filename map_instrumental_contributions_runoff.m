% figure set up
clr = cbrewer('div','RdBu',14); clr(clr<0) = 0; clr(clr>1) = 1;

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 0 4.75 9];

%% MsTMIP early instrumental
load ./data/MsTMIP_WaterBudget;
load ./data/MRB_subregions_GP;
[~, ny, nx, nk] = size(ET_RG1);
states = shaperead('usastatehi','UseGeoCoords', true);
latlim=[36 50];
lonlim=[-115 -89];

% convert from meters per year to mm per year
P = P * 1000;
ET_SG3 = ET_SG3 * 1000;
ET_SG2 = ET_SG2 * 1000;
ET_SG1 = ET_SG1 * 1000;
ET_RG1 = ET_RG1 * 1000;

% get subbasin number
[LON, LAT] = meshgrid(lon, lat);
LatLon = [reshape(LAT, [], 1) reshape(LON, [], 1)];
idx = zeros(size(LatLon,1),1);
[IN, ON] = inpolygon(LatLon(:,1), LatLon(:,2), SR1(1).Lat, SR1(1).Lon);
idx(IN | ON) = 1;
[IN, ON] = inpolygon(LatLon(:,1), LatLon(:,2), SR2(2).Lat, SR2(2).Lon);
idx(IN | ON) = 2;
[IN, ON] = inpolygon(LatLon(:,1), LatLon(:,2), SR3(2).Lat, SR3(2).Lon);
idx(IN | ON) = 3;
[IN, ON] = inpolygon(LatLon(:,1), LatLon(:,2), SR4(2).Lat, SR4(2).Lon);
idx(IN | ON) = 4;
[IN, ON] = inpolygon(LatLon(:,1), LatLon(:,2), SR5(2).Lat, SR5(2).Lon);
idx(IN | ON) = 5;
[IN, ON] = inpolygon(LatLon(:,1), LatLon(:,2), SR6(2).Lat, SR6(2).Lon);
idx(IN | ON) = 6;

MRBidx = reshape(idx, ny, nx);
clear LAT LON LatLon IN ON idx e;

% get contributions
idx = year>=1981 & year<=2010;
nt = length(year); 
Pmean = repmat(mean(P(1:30,:,:),'omitnan'),nt,1,1,nk);

WB_RG1 = Pmean - ET_RG1;
WB_SG1 = repmat(P, 1, 1, 1, nk) - ET_SG1;
WB_SG2 = repmat(P, 1, 1, 1, nk) - ET_SG2;
WB_SG3 = repmat(P, 1, 1, 1, nk) - ET_SG3;

Clim = WB_SG1 - WB_RG1;
LULCC = WB_SG2 - WB_SG1;
CO2 = WB_SG3 - WB_SG2;
All = WB_SG3 - WB_RG1;

Clim = squeeze(mean(Clim(idx,:,:,:),[1 4]) - mean(Clim(year>=1931 & year<=1960,:,:,:),[1 4])); Clim(MRBidx==0) = NaN;
LULCC = squeeze(mean(LULCC(idx,:,:,:),[1 4]) - mean(LULCC(year>=1931 & year<=1960,:,:,:),[1 4])); LULCC(MRBidx==0) = NaN;
CO2 = squeeze(mean(CO2(idx,:,:,:),[1 4]) - mean(CO2(year>=1931 & year<=1960,:,:,:),[1 4])); CO2(MRBidx==0) = NaN;
All = squeeze(mean(All(idx,:,:,:),[1 4]) - mean(All(year>=1931 & year<=1960,:,:,:),[1 4])); All(MRBidx==0) = NaN;

% log2
Clim(Clim>0) = log2(Clim(Clim>0)); Clim(Clim<0) = -log2(abs(Clim(Clim<0)));
LULCC(LULCC>0) = log2(LULCC(LULCC>0)); LULCC(LULCC<0) = -log2(abs(LULCC(LULCC<0)));
CO2(CO2>0) = log2(CO2(CO2>0)); CO2(CO2<0) = -log2(abs(CO2(CO2<0)));

% Map
subplot(11,4,[1 2 5 6])
axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'off','PLineLocation',4,'MLineLocation',6,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.3,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','MLabelParallel','north',...
        'FLineWidth',1, 'GColor',[0.5 0.5 0.5], 'FontSize',8)
axis off;
axis image;
surfm(lat-0.25, lon-0.25, CO2);
caxis([-7 7])
colormap(gca, clr)
geoshow(states,'FaceColor','none','EdgeColor',[0.5 0.5 0.5],'LineWidth',0.5)
geoshow(SR1(1,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR2(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR3(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR4(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR5(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR6(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
ax = gca;
subplotsqueeze(ax, 1.1)
ax.Position(1) = 0.05;
ax.Position(2) = 0.82;
text(-0.18, 0.83, 'a', 'FontSize',12)
ttl = title('\DeltaR_{s,CO2}','FontSize',9,'FontWeight','normal');
ttl.Position(2) = 0.83;

subplot(11,4,[3 4 7 8])
axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'off','PLineLocation',4,'MLineLocation',6,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.3,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','MLabelParallel','north',...
        'FLineWidth',1, 'GColor',[0.5 0.5 0.5], 'FontSize',8)
axis off;
axis image;
surfm(lat-0.25, lon-0.25, LULCC);
caxis([-7 7])
colormap(gca, clr)
geoshow(states,'FaceColor','none','EdgeColor',[0.5 0.5 0.5],'LineWidth',0.5)
geoshow(SR1(1,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR2(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR3(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR4(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR5(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR6(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
ax = gca;
subplotsqueeze(ax, 1.1)
ax.Position(1) = 0.48;
ax.Position(2) = 0.82;
text(-0.18, 0.83, 'b', 'FontSize',12)
ttl = title('\DeltaR_{s,LULCC}','FontSize',9,'FontWeight','normal');
ttl.Position(2) = 0.83;

subplot(11,4,[9 10 13 14])
axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'off','PLineLocation',4,'MLineLocation',6,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.3,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','MLabelParallel','north',...
        'FLineWidth',1, 'GColor',[0.5 0.5 0.5], 'FontSize',8)
axis off;
axis image;
surfm(lat-0.25, lon-0.25, Clim);
caxis([-7 7])
colormap(gca, clr)
geoshow(states,'FaceColor','none','EdgeColor',[0.5 0.5 0.5],'LineWidth',0.5)
geoshow(SR1(1,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR2(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR3(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR4(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR5(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR6(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
ax = gca;
subplotsqueeze(ax, 1.1)
ax.Position(1) = 0.06;
ax.Position(2) = 0.64;
text(-0.18, 0.83, 'c', 'FontSize',12)
ttl = title('\DeltaR_{s,climate} (MsTMIP)','FontSize',9,'FontWeight','normal');
ttl.Position(2) = 0.83;


%% map climate contributions runoff

load ./data/McCabeWilliams_WaterBudget;
syear = 1981;
eyear = 2010;
[~, ny, nx] = size(Rs_S3_WY);
cellsize = 0.25;
days_in_month = [31 28 31 30 31 30 31 31 30 31 30 31];
windowSize = 12;
bx = ones(1,windowSize);
a = 1;


%% get subbasin number
[LON, LAT] = meshgrid(lon, lat);
LatLon = [reshape(LAT, [], 1) reshape(LON, [], 1)];
idx = zeros(size(LatLon,1),1);
[IN, ON] = inpolygon(LatLon(:,1), LatLon(:,2), SR1(1).Lat, SR1(1).Lon);
idx(IN | ON) = 1;
[IN, ON] = inpolygon(LatLon(:,1), LatLon(:,2), SR2(2).Lat, SR2(2).Lon);
idx(IN | ON) = 2;
[IN, ON] = inpolygon(LatLon(:,1), LatLon(:,2), SR3(2).Lat, SR3(2).Lon);
idx(IN | ON) = 3;
[IN, ON] = inpolygon(LatLon(:,1), LatLon(:,2), SR4(2).Lat, SR4(2).Lon);
idx(IN | ON) = 4;
[IN, ON] = inpolygon(LatLon(:,1), LatLon(:,2), SR5(2).Lat, SR5(2).Lon);
idx(IN | ON) = 5;
[IN, ON] = inpolygon(LatLon(:,1), LatLon(:,2), SR6(2).Lat, SR6(2).Lon);
idx(IN | ON) = 6;

MRBidx = reshape(idx, ny, nx);

[x1, y1] = centroid(polyshape(SR1(1).Lon, SR1(1).Lat));
[x2, y2] = centroid(polyshape(SR2(2).Lon, SR2(2).Lat));
[x3, y3] = centroid(polyshape(SR3(2).Lon, SR3(2).Lat));
[x4, y4] = centroid(polyshape(SR4(2).Lon, SR4(2).Lat));
[x5, y5] = centroid(polyshape(SR5(2).Lon, SR5(2).Lat));
[x6, y6] = centroid(polyshape(SR6(2).Lon, SR6(2).Lat));

e = referenceEllipsoid('World Geodetic System 1984');
area = areaquad(reshape(LAT-(cellsize/2),[],1),reshape(LON-(cellsize/2),[],1),reshape(LAT+(cellsize/2),[],1),reshape(LON+(cellsize/2),[],1),e);
area = reshape(area, ny, nx); 

clear LAT LON LatLon IN ON idx e;

%% get monthly contributions
yr = reshape(repmat(year, 12, 1), [], 1);
mos = repmat(1:12, 1, length(year))';
nt = length(yr); 
clat = readmatrix('./modeling/lat.txt');
clon = readmatrix('./modeling/lon.txt');
s0 = readmatrix("./modeling/wb.lmrb.awc.runoff.hamonpet.s0", 'FileType','text');
s1 = readmatrix("./modeling/wb.lmrb.awc.runoff.hamonpet.s1", 'FileType','text');
s2 = readmatrix("./modeling/wb.lmrb.awc.runoff.hamonpet.s2", 'FileType','text');
s3 = readmatrix("./modeling/wb.lmrb.awc.runoff.hamonpet.s3", 'FileType','text');

Rs_S0 = NaN(nt, ny, nx);
Rs_S1 = NaN(nt, ny, nx);
Rs_S2 = NaN(nt, ny, nx);
Rs_S3 = NaN(nt, ny, nx);

for i = 1:length(clat)
    xi = find(lon == clon(i));
    yi = find(lat == clat(i));

    Rs_S0(yr >= s0(1,1), yi, xi) = s0(:,i+2);
    Rs_S1(yr >= s1(1,1), yi, xi) = s1(:,i+2);
    Rs_S2(yr >= s2(1,1), yi, xi) = s2(:,i+2);
    Rs_S3(yr >= s3(1,1), yi, xi) = s3(:,i+2);
end

%% map contributions
% Natural vs. Anthropogenic totals
Qnat = squeeze(mean(Rs_S0_WY(year>=syear & year<=eyear, :, :) - Rs_S0_WY(year>=1931 & year<=1960, :, :)));
Qanth = squeeze(mean(Rs_S3_WY(year>=syear & year<=eyear, :, :) - Rs_S3_WY(year>=1931 & year<=1960, :, :))) - Qnat;
Qall = Qnat + Qanth;

% log2
Qnat(Qnat>0) = log2(Qnat(Qnat>0)); Qnat(Qnat<0) = -log2(abs(Qnat(Qnat<0)));
Qanth(Qanth>0) = log2(Qanth(Qanth>0)); Qanth(Qanth<0) = -log2(abs(Qanth(Qanth<0)));
Qall(Qall>0) = log2(Qall(Qall>0)); Qall(Qall<0) = -log2(abs(Qall(Qall<0)));


subplot(11,4,[11 12 15 16])
axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'off','PLineLocation',4,'MLineLocation',6,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.3,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','MLabelParallel','north',...
        'FLineWidth',1, 'GColor',[0.5 0.5 0.5], 'FontSize',8)
axis off;
axis image;
surfm(lat, lon, Qnat+Qanth);
caxis([-7 7])
colormap(gca, clr)
geoshow(states,'FaceColor','none','EdgeColor',[0.5 0.5 0.5],'LineWidth',0.5)
geoshow(SR1(1,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR2(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR3(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR4(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR5(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR6(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
ax = gca;
subplotsqueeze(ax, 1.1)
ax.Position(1) = 0.48;
ax.Position(2) = 0.64;
text(-0.18, 0.83, 'd', 'FontSize',12)
ttl = title('\DeltaR_{s,climate} (MW11)','FontSize',9,'FontWeight','normal');
ttl.Position(2) = 0.83;
annotation("line",[0.05 0.89],[0.625 0.625], 'LineWidth',1.5)
annotation("line",[0.05 0.05],[0.625 0.61], 'LineWidth',1.5)
annotation("line",[0.89 0.89],[0.625 0.61], 'LineWidth',1.5)
annotation("line",[0.68 0.68],[0.625 0.637], 'LineWidth',1.5)

subplot(11,4,[17 18 21 22])
axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'off','PLineLocation',4,'MLineLocation',6,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.3,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','MLabelParallel','north',...
        'FLineWidth',1, 'GColor',[0.5 0.5 0.5], 'FontSize',8)
axis off;
axis image;
surfm(lat, lon, Qnat);
caxis([-7 7])
colormap(gca, clr)
geoshow(states,'FaceColor','none','EdgeColor',[0.5 0.5 0.5],'LineWidth',0.5)
geoshow(SR1(1,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR2(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR3(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR4(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR5(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR6(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
ax = gca;
subplotsqueeze(ax, 1.1)
ax.Position(1) = 0.06;
ax.Position(2) = 0.46;
text(-0.18, 0.83, 'e', 'FontSize',12)
ttl = title('\DeltaR_{s,climate} (natural)','FontSize',9,'FontWeight','normal');
ttl.Position(2) = 0.83;

subplot(11,4,[19 20 23 24])
axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'off','PLineLocation',4,'MLineLocation',6,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.3,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','MLabelParallel','north',...
        'FLineWidth',1, 'GColor',[0.5 0.5 0.5], 'FontSize',8)
axis off;
axis image;
surfm(lat, lon, Qanth);
caxis([-7 7])
colormap(gca, clr)
geoshow(states,'FaceColor','none','EdgeColor',[0.5 0.5 0.5],'LineWidth',0.5)
geoshow(SR1(1,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR2(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR3(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR4(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR5(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR6(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
ax = gca;
subplotsqueeze(ax, 1.1)
ax.Position(1) = 0.48;
ax.Position(2) = 0.46;
text(-0.18, 0.83, 'f', 'FontSize',12)
ttl = title('\DeltaR_{s,climate} (anthropogenic)','FontSize',9,'FontWeight','normal');
ttl.Position(2) = 0.83;
annotation("line",[0.03 0.95],[0.44 0.44], 'LineWidth',1.5)
annotation("line",[0.03 0.03],[0.44 0.425], 'LineWidth',1.5)
annotation("line",[0.95 0.95],[0.44 0.425], 'LineWidth',1.5)
annotation("line",[0.68 0.68],[0.44 0.455], 'LineWidth',1.5)

% subcomponents
h1 = axes('Parent', gcf, 'Position', [0.03 0.28 0.3 0.15]);
set(h1, 'Color','w')
axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'off','PLineLocation',4,'MLineLocation',6,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.3,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','MLabelParallel','north',...
        'FLineWidth',1, 'GColor',[0.5 0.5 0.5], 'FontSize',8)
axis off;
axis image;
temp = squeeze(mean(Rs_S1_WY(year>=syear & year<=eyear, :, :) - Rs_S0_WY(year>=syear & year<=eyear, :, :)));
temp(temp>0) = log2(temp(temp>0)); temp(temp<0) = -log2(abs(temp(temp<0)));
surfm(lat, lon, temp);
caxis([-7 7])
colormap(gca, clr)
geoshow(states,'FaceColor','none','EdgeColor',[0.5 0.5 0.5],'LineWidth',0.5)
geoshow(SR1(1,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR2(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR3(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR4(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR5(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR6(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
textm(y1, x1, '1', 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontWeight','bold', 'FontSize',12)
textm(y2, x2, '2', 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontWeight','bold', 'FontSize',12)
textm(y3, x3, '3', 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontWeight','bold', 'FontSize',12)
textm(y4, x4, '4', 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontWeight','bold', 'FontSize',12)
textm(y5, x5, '5', 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontWeight','bold', 'FontSize',12)
textm(y6, x6, '6', 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontWeight','bold', 'FontSize',12)
text(-0.18, 0.84, 'g', 'FontSize',12)
ttl = title('\DeltaR_{s,climate} (T_{avg})','FontSize',9,'FontWeight','normal');
ttl.Position(2) = 0.83;

h1 = axes('Parent', gcf, 'Position', [0.35 0.28 0.3 0.15]);
set(h1, 'Color','w')
axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'off','PLineLocation',4,'MLineLocation',6,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.3,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','MLabelParallel','north',...
        'FLineWidth',1, 'GColor',[0.5 0.5 0.5], 'FontSize',8)
axis off;
axis image;
temp = squeeze(mean(Rs_S3_WY(year>=syear & year<=eyear, :, :) - Rs_S2_WY(year>=syear & year<=eyear, :, :)));
temp(temp>0) = log2(temp(temp>0)); temp(temp<0) = -log2(abs(temp(temp<0)));
surfm(lat, lon, temp);
caxis([-7 7])
colormap(gca, clr)
geoshow(states,'FaceColor','none','EdgeColor',[0.5 0.5 0.5],'LineWidth',0.5)
geoshow(SR1(1,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR2(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR3(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR4(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR5(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR6(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
text(-0.18, 0.84, 'h', 'FontSize',12)
ttl = title('\DeltaR_{s,climate} (P)','FontSize',9,'FontWeight','normal');
ttl.Position(2) = 0.83;

h1 = axes('Parent', gcf, 'Position', [0.67 0.28 0.3 0.15]);
set(h1, 'Color','w')
axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'off','PLineLocation',4,'MLineLocation',6,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.3,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','MLabelParallel','north',...
        'FLineWidth',1, 'GColor',[0.5 0.5 0.5], 'FontSize',8)
axis off;
axis image;
temp = squeeze(mean(Rs_S2_WY(year>=syear & year<=eyear, :, :) - Rs_S1_WY(year>=syear & year<=eyear, :, :)));
surfm(lat, lon, temp);
caxis([-7 7])
colormap(gca, clr)
geoshow(states,'FaceColor','none','EdgeColor',[0.5 0.5 0.5],'LineWidth',0.5)
geoshow(SR1(1,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR2(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR3(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR4(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR5(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR6(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
text(-0.18, 0.84, 'i', 'FontSize',12)
ttl = title('\DeltaR_{s,climate} (PET)','FontSize',9,'FontWeight','normal');
ttl.Position(2) = 0.83;

cb = colorbar('southoutside');
cb.Position = [0.03 0.28 0.94 0.01];
cb.Ticks = -7:1:7;
cb.TickLength = 0.018;
cb.FontSize = 7;
cb.TickLabels = {'-128','-64','-32','-16','-8','-4','-2','0','2','4','8','16','32','64','128'};
ylabel(cb, '\DeltaR_{s} (mm)', 'FontSize',8)

%% Contributions by season / subbasin
clr = wesanderson('fantasticfox1');
clr = clr(1:3,:);

subplot(11,4,37:44)
Q_S0 = reshape(sum((0.001 .* Rs_S0(:, MRBidx > 0) .* repmat(area(MRBidx > 0)', size(Rs_S0, 1), 1)), 2), 12, []) ./ (repmat(days_in_month', 1, length(year)) * 24 * 60 * 60);
Q_S1 = reshape(sum((0.001 .* Rs_S1(:, MRBidx > 0) .* repmat(area(MRBidx > 0)', size(Rs_S1, 1), 1)), 2), 12, []) ./ (repmat(days_in_month', 1, length(year)) * 24 * 60 * 60);
Q_S2 = reshape(sum((0.001 .* Rs_S2(:, MRBidx > 0) .* repmat(area(MRBidx > 0)', size(Rs_S2, 1), 1)), 2), 12, []) ./ (repmat(days_in_month', 1, length(year)) * 24 * 60 * 60);
Q_S3 = reshape(sum((0.001 .* Rs_S3(:, MRBidx > 0) .* repmat(area(MRBidx > 0)', size(Rs_S3, 1), 1)), 2), 12, []) ./ (repmat(days_in_month', 1, length(year)) * 24 * 60 * 60);
dQ_tavg = Q_S1 - Q_S0;
dQ_pet = Q_S2 - Q_S1;
dQ_p = Q_S3 - Q_S2;
dat = [mean(dQ_tavg(:,year>=syear & year<=eyear), 2) - mean(dQ_tavg(:,year>=1931 & year<=1960), 2) ...
    mean(dQ_pet(:,year>=syear & year<=eyear), 2) - mean(dQ_pet(:,year>=1931 & year<=1960), 2) ...
    mean(dQ_p(:,year>=syear & year<=eyear), 2) - mean(dQ_p(:,year>=1931 & year<=1960), 2)];
Q_S0 = filter(bx, a, sum((0.001 .* Rs_S0(:, MRBidx > 0) .* repmat(area(MRBidx > 0)', size(Rs_S0, 1), 1)), 2), [], 1); Q_S0 = Q_S0(mos==9) / (365 * 24 * 60 * 60);
Q_S1 = filter(bx, a, sum((0.001 .* Rs_S1(:, MRBidx > 0) .* repmat(area(MRBidx > 0)', size(Rs_S1, 1), 1)), 2), [], 1); Q_S1 = Q_S1(mos==9) / (365 * 24 * 60 * 60);
Q_S2 = filter(bx, a, sum((0.001 .* Rs_S2(:, MRBidx > 0) .* repmat(area(MRBidx > 0)', size(Rs_S2, 1), 1)), 2), [], 1); Q_S2 = Q_S2(mos==9) / (365 * 24 * 60 * 60);
Q_S3 = filter(bx, a, sum((0.001 .* Rs_S3(:, MRBidx > 0) .* repmat(area(MRBidx > 0)', size(Rs_S3, 1), 1)), 2), [], 1); Q_S3 = Q_S3(mos==9) / (365 * 24 * 60 * 60);
dQ_tavg = Q_S1 - Q_S0;
dQ_pet = Q_S2 - Q_S1;
dQ_p = Q_S3 - Q_S2;
dat(13,:) = [mean(dQ_tavg(year>=syear & year<=eyear)) - mean(dQ_tavg(year>=1931 & year<=1960)) ...
    mean(dQ_pet(year>=syear & year<=eyear)) - mean(dQ_pet(year>=1931 & year<=1960)) ...
    mean(dQ_p(year>=syear & year<=eyear)) - mean(dQ_p(year>=1931 & year<=1960))];
b = bar(1:13,dat,'stacked');
b(1).FaceColor = clr(1,:);
b(2).FaceColor = clr(2,:);
b(3).FaceColor = clr(3,:);
ax = gca;
set(ax, 'XLim',[0.5 13.5], 'TickDir','out', 'XTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D', 'WY'})
box off;
hold on;
scatter(1:13, sum(dat,2), 20, "black","filled")
ax.Position(1) = 0.12;
ax.Position(2) = 0.1;
ax.Position(3) = 0.84;
ylb = ylabel(ax, '\DeltaQ (m^{3} s^{-1})', 'FontSize',8);
ylb.Position(1) = -0.7;
ylim = get(ax, 'YLim');
lgd = legend('T_{avg}','PET','P', 'Location','northoutside', 'Orientation', 'horizontal');
lgd.Position(2) = 0.21;
lgd.Position(1) = 0.5;
legend('boxoff')
text(-1.1,ylim(2),'j','FontSize',12)

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300',['./output/map-instrumental-contributions-runoff-',num2str(syear),'-',num2str(eyear),'.tif'])
close all;




