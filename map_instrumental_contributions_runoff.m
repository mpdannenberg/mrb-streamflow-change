% figure set up
clr = cbrewer('div','RdBu',14); clr(clr<0) = 0; clr(clr>1) = 1;

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 0 4.75 6.75];

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
subplot(9,4,[1 2 5 6])
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
ax.Position(2) = 0.8;
text(-0.18, 0.83, 'A', 'FontSize',12, 'FontWeight','bold')
ttl = title('\DeltaR_{s,CO2}','FontSize',9,'FontWeight','normal');
ttl.Position(2) = 0.83;

subplot(9,4,[3 4 7 8])
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
ax.Position(2) = 0.8;
text(-0.18, 0.83, 'B', 'FontSize',12, 'FontWeight','bold')
ttl = title('\DeltaR_{s,LULCC}','FontSize',9,'FontWeight','normal');
ttl.Position(2) = 0.83;

subplot(9,4,[9 10 13 14])
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
ax.Position(2) = 0.58;
text(-0.18, 0.83, 'C', 'FontSize',12, 'FontWeight','bold')
ttl = title('\DeltaR_{s,climate} (MsTMIP)','FontSize',9,'FontWeight','normal');
ttl.Position(2) = 0.83;


%% map climate contributions runoff
load ./data/McCabeWilliams_WaterBudget;
syear = 1981;
eyear = 2010;
[~, ny, nx] = size(Rs_S3_WY);

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
clear LAT LON LatLon IN ON idx;

%% map contributions
% Natural vs. Anthropogenic totals
Qnat = squeeze(mean(Rs_S0_WY(year>=syear & year<=eyear, :, :) - Rs_S0_WY(year>=1931 & year<=1960, :, :)));
Qanth = squeeze(mean(Rs_S3_WY(year>=syear & year<=eyear, :, :) - Rs_S3_WY(year>=1931 & year<=1960, :, :))) - Qnat;
Qall = Qnat + Qanth;

% log2
Qnat(Qnat>0) = log2(Qnat(Qnat>0)); Qnat(Qnat<0) = -log2(abs(Qnat(Qnat<0)));
Qanth(Qanth>0) = log2(Qanth(Qanth>0)); Qanth(Qanth<0) = -log2(abs(Qanth(Qanth<0)));
Qall(Qall>0) = log2(Qall(Qall>0)); Qall(Qall<0) = -log2(abs(Qall(Qall<0)));

subplot(9,4,[11 12 15 16])
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
ax.Position(2) = 0.58;
text(-0.18, 0.83, 'D', 'FontSize',12, 'FontWeight','bold')
ttl = title('\DeltaR_{s,climate} (MW11)','FontSize',9,'FontWeight','normal');
ttl.Position(2) = 0.83;
annotation("line",[0.05 0.89],[0.625-0.06 0.625-0.06], 'LineWidth',1.5)
annotation("line",[0.05 0.05],[0.625-0.06 0.61-0.06], 'LineWidth',1.5)
annotation("line",[0.89 0.89],[0.625-0.06 0.61-0.06], 'LineWidth',1.5)
annotation("line",[0.68 0.68],[0.625-0.06 0.637-0.06], 'LineWidth',1.5)

subplot(9,4,[17 18 21 22])
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
ax.Position(2) = 0.36;
text(-0.18, 0.83, 'E', 'FontSize',12, 'FontWeight','bold')
ttl = title('\DeltaR_{s,climate} (natural)','FontSize',9,'FontWeight','normal');
ttl.Position(2) = 0.83;

subplot(9,4,[19 20 23 24])
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
ax.Position(2) = 0.36;
text(-0.18, 0.83, 'F', 'FontSize',12, 'FontWeight','bold')
ttl = title('\DeltaR_{s,climate} (anthropogenic)','FontSize',9,'FontWeight','normal');
ttl.Position(2) = 0.83;
annotation("line",[0.03 0.95],[0.44-0.1 0.44-0.1], 'LineWidth',1.5)
annotation("line",[0.03 0.03],[0.44-0.1 0.425-0.1], 'LineWidth',1.5)
annotation("line",[0.95 0.95],[0.44-0.1 0.425-0.1], 'LineWidth',1.5)
annotation("line",[0.68 0.68],[0.44-0.1 0.455-0.1], 'LineWidth',1.5)

% subcomponents
h1 = axes('Parent', gcf, 'Position', [0.03 0.15 0.3 0.15]);
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
text(-0.18, 0.84, 'G', 'FontSize',12, 'FontWeight','bold')
ttl = title('\DeltaR_{s,climate} (T_{avg})','FontSize',9,'FontWeight','normal');
ttl.Position(2) = 0.83;

h1 = axes('Parent', gcf, 'Position', [0.35 0.15 0.3 0.15]);
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
text(-0.18, 0.84, 'H', 'FontSize',12, 'FontWeight','bold')
ttl = title('\DeltaR_{s,climate} (P)','FontSize',9,'FontWeight','normal');
ttl.Position(2) = 0.83;

h1 = axes('Parent', gcf, 'Position', [0.67 0.15 0.3 0.15]);
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
text(-0.18, 0.84, 'I', 'FontSize',12, 'FontWeight','bold')
ttl = title('\DeltaR_{s,climate} (PET)','FontSize',9,'FontWeight','normal');
ttl.Position(2) = 0.83;

cb = colorbar('southoutside');
cb.Position = [0.03 0.1 0.94 0.02];
cb.Ticks = -7:1:7;
cb.TickLength = 0.03;
cb.FontSize = 8;
cb.TickLabels = {'-128','-64','-32','-16','-8','-4','-2','0','2','4','8','16','32','64','128'};
ylabel(cb, '\DeltaR_{s} (mm)', 'FontSize',10)

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300',['./output/map-instrumental-contributions-runoff-',num2str(syear),'-',num2str(eyear),'.tif'])
close all;

