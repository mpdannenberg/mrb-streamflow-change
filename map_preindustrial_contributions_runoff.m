% map climate contributions runoff

load ./data/MsTMIP_WaterBudget;
load ./data/MRB_subregions_GP;
syear = 1901;
eyear = 2010;
[~, ny, nx, nk] = size(ET_RG1);
states = shaperead('usastatehi','UseGeoCoords', true);
latlim=[36 50];
lonlim=[-115 -89];

%% convert from meters per year to mm per year
P = P * 1000;
ET_SG3 = ET_SG3 * 1000;
ET_SG2 = ET_SG2 * 1000;
ET_SG1 = ET_SG1 * 1000;
ET_RG1 = ET_RG1 * 1000;

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
area = areaquad(reshape(LAT-0.25,[],1),reshape(LON-0.25,[],1),reshape(LAT+0.25,[],1),reshape(LON+0.25,[],1),e);
area = reshape(area, ny, nx); 

clear LAT LON LatLon IN ON idx e;

%% get contributions
idx = year>=1981 & year<=2010;
nt = length(year); 
Pmean = repmat(mean(P(1:30,:,:),'omitnan'),nt,1,1,nk);

WB_RG1 = Pmean - ET_RG1;
WB_SG1 = repmat(P, 1, 1, 1, nk) - ET_SG1;
WB_SG2 = repmat(P, 1, 1, 1, nk) - ET_SG2;
WB_SG3 = repmat(P, 1, 1, 1, nk) - ET_SG3;

Clim = squeeze(mean(mean(WB_SG1(idx,:,:,:) - WB_RG1(idx,:,:,:), 4), 1)); Clim(MRBidx==0) = NaN;
LULCC = squeeze(mean(mean(WB_SG2(idx,:,:,:) - WB_SG1(idx,:,:,:), 4), 1)); LULCC(MRBidx==0) = NaN;
CO2 = squeeze(mean(mean(WB_SG3(idx,:,:,:) - WB_SG2(idx,:,:,:), 4), 1)); CO2(MRBidx==0) = NaN;
All = squeeze(mean(mean(WB_SG3(idx,:,:,:) - WB_RG1(idx,:,:,:), 4), 1)); All(MRBidx==0) = NaN;

% log2
Clim(Clim>0) = log2(Clim(Clim>0)); Clim(Clim<0) = -log2(abs(Clim(Clim<0)));
LULCC(LULCC>0) = log2(LULCC(LULCC>0)); LULCC(LULCC<0) = -log2(abs(LULCC(LULCC<0)));
CO2(CO2>0) = log2(CO2(CO2>0)); CO2(CO2<0) = -log2(abs(CO2(CO2<0)));
All(All>0) = log2(All(All>0)); All(All<0) = -log2(abs(All(All<0)));

%% map contributions
clr = cbrewer('div','RdBu',14); clr(clr<0) = 0; clr(clr>1) = 1;
[LON, LAT] = meshgrid(lon,lat);

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6.5 7.3];

subplot(7,2,[1 3])
axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'off','PLineLocation',4,'MLineLocation',6,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.3,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','MLabelParallel','north',...
        'FLineWidth',1, 'GColor',[0.5 0.5 0.5], 'FontSize',8)
axis off;
axis image;
surfm(LAT-0.25, LON-0.25, All);
caxis([-7 7])
colormap(gca, clr)
geoshow(states,'FaceColor','none','EdgeColor',[0.5 0.5 0.5],'LineWidth',0.5)
geoshow(SR1(1,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR2(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR3(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR4(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR5(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR6(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
ttl = title('\DeltaR_{s}', 'FontSize',9,'FontWeight','normal'); ttl.Position(2) = 0.84;
ax = gca;
ax.Position(1) = 0.06;
ax.Position(2) = 0.76;
text(-0.16,0.83,'a','FontSize',12)

subplot(7,2,[2 4])
axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'off','PLineLocation',4,'MLineLocation',6,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.3,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','MLabelParallel','north',...
        'FLineWidth',1, 'GColor',[0.5 0.5 0.5], 'FontSize',8)
axis off;
axis image;
surfm(LAT-0.25, LON-0.25, Clim);
caxis([-7 7])
colormap(gca, clr)
geoshow(states,'FaceColor','none','EdgeColor',[0.5 0.5 0.5],'LineWidth',0.5)
geoshow(SR1(1,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR2(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR3(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR4(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR5(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR6(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
ttl = title('\DeltaR_{s,climate}', 'FontSize',9,'FontWeight','normal'); ttl.Position(2) = 0.84;
ax = gca;
ax.Position(1) = 0.48;
ax.Position(2) = 0.76;
text(-0.16,0.83,'b','FontSize',12)

subplot(7,2,[5 7])
axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'off','PLineLocation',4,'MLineLocation',6,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.3,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','MLabelParallel','north',...
        'FLineWidth',1, 'GColor',[0.5 0.5 0.5], 'FontSize',8)
axis off;
axis image;
surfm(LAT-0.25, LON-0.25, LULCC);
caxis([-7 7])
colormap(gca, clr)
geoshow(states,'FaceColor','none','EdgeColor',[0.5 0.5 0.5],'LineWidth',0.5)
geoshow(SR1(1,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR2(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR3(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR4(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR5(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR6(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
ttl = title('\DeltaR_{s,LULCC}', 'FontSize',9, 'FontWeight','normal'); ttl.Position(2) = 0.84;
ax = gca;
ax.Position(1) = 0.06;
ax.Position(2) = 0.5;
text(-0.16,0.83,'c','FontSize',12)

subplot(7,2,[6 8])
axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'off','PLineLocation',4,'MLineLocation',6,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.3,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','MLabelParallel','north',...
        'FLineWidth',1, 'GColor',[0.5 0.5 0.5], 'FontSize',8)
axis off;
axis image;
surfm(LAT-0.25, LON-0.25, CO2);
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
ttl = title('\DeltaR_{s,CO2}', 'FontSize',9, 'FontWeight','normal'); ttl.Position(2) = 0.84;
ax = gca;
ax.Position(1) = 0.48;
ax.Position(2) = 0.5;
text(-0.16,0.83,'d','FontSize',12)

cb = colorbar('eastoutside');
cb.Position = [0.85 0.52 0.03 0.44];
cb.Ticks = -7:1:7;
cb.TickLabels = {'-128','-64','-32','-16','-8','-4','-2','0','2','4','8','16','32','64','128'};
cb.TickLength = 0.062;
ylabel(cb, '\DeltaR_{s} (mm)', 'FontSize',10)

%% make contribution plots, by region, below that one
clr = wesanderson('fantasticfox1');
clr = clr([2 1 3],:);
Q_RG1 = NaN(nt, 7, nk);
Q_SG1 = NaN(nt, 7, nk);
Q_SG2 = NaN(nt, 7, nk);
Q_SG3 = NaN(nt, 7, nk);
for i = 1:nk
    Rs_RG1 = WB_RG1(:,:,:,i);
    Rs_SG1 = WB_SG1(:,:,:,i);
    Rs_SG2 = WB_SG2(:,:,:,i);
    Rs_SG3 = WB_SG3(:,:,:,i);

    % All regions
    Q_RG1(:,1,i) = sum((0.001 .* Rs_RG1(:, MRBidx > 0) .* repmat(area(MRBidx > 0)', size(Rs_RG1, 1), 1)) / (365 * 24 * 60 * 60), 2);
    Q_SG1(:,1,i) = sum((0.001 .* Rs_SG1(:, MRBidx > 0) .* repmat(area(MRBidx > 0)', size(Rs_SG1, 1), 1)) / (365 * 24 * 60 * 60), 2);
    Q_SG2(:,1,i) = sum((0.001 .* Rs_SG2(:, MRBidx > 0) .* repmat(area(MRBidx > 0)', size(Rs_SG2, 1), 1)) / (365 * 24 * 60 * 60), 2);
    Q_SG3(:,1,i) = sum((0.001 .* Rs_SG3(:, MRBidx > 0) .* repmat(area(MRBidx > 0)', size(Rs_SG3, 1), 1)) / (365 * 24 * 60 * 60), 2);

    % Region 1
    Q_RG1(:,2,i) = sum((0.001 .* Rs_RG1(:, MRBidx == 1) .* repmat(area(MRBidx == 1)', size(Rs_RG1, 1), 1)) / (365 * 24 * 60 * 60), 2);
    Q_SG1(:,2,i) = sum((0.001 .* Rs_SG1(:, MRBidx == 1) .* repmat(area(MRBidx == 1)', size(Rs_SG1, 1), 1)) / (365 * 24 * 60 * 60), 2);
    Q_SG2(:,2,i) = sum((0.001 .* Rs_SG2(:, MRBidx == 1) .* repmat(area(MRBidx == 1)', size(Rs_SG2, 1), 1)) / (365 * 24 * 60 * 60), 2);
    Q_SG3(:,2,i) = sum((0.001 .* Rs_SG3(:, MRBidx == 1) .* repmat(area(MRBidx == 1)', size(Rs_SG3, 1), 1)) / (365 * 24 * 60 * 60), 2);

    % Region 2
    Q_RG1(:,3,i) = sum((0.001 .* Rs_RG1(:, MRBidx == 2) .* repmat(area(MRBidx == 2)', size(Rs_RG1, 1), 1)) / (365 * 24 * 60 * 60), 2);
    Q_SG1(:,3,i) = sum((0.001 .* Rs_SG1(:, MRBidx == 2) .* repmat(area(MRBidx == 2)', size(Rs_SG1, 1), 1)) / (365 * 24 * 60 * 60), 2);
    Q_SG2(:,3,i) = sum((0.001 .* Rs_SG2(:, MRBidx == 2) .* repmat(area(MRBidx == 2)', size(Rs_SG2, 1), 1)) / (365 * 24 * 60 * 60), 2);
    Q_SG3(:,3,i) = sum((0.001 .* Rs_SG3(:, MRBidx == 2) .* repmat(area(MRBidx == 2)', size(Rs_SG3, 1), 1)) / (365 * 24 * 60 * 60), 2);

    % Region 3
    Q_RG1(:,4,i) = sum((0.001 .* Rs_RG1(:, MRBidx == 3) .* repmat(area(MRBidx == 3)', size(Rs_RG1, 1), 1)) / (365 * 24 * 60 * 60), 2);
    Q_SG1(:,4,i) = sum((0.001 .* Rs_SG1(:, MRBidx == 3) .* repmat(area(MRBidx == 3)', size(Rs_SG1, 1), 1)) / (365 * 24 * 60 * 60), 2);
    Q_SG2(:,4,i) = sum((0.001 .* Rs_SG2(:, MRBidx == 3) .* repmat(area(MRBidx == 3)', size(Rs_SG2, 1), 1)) / (365 * 24 * 60 * 60), 2);
    Q_SG3(:,4,i) = sum((0.001 .* Rs_SG3(:, MRBidx == 3) .* repmat(area(MRBidx == 3)', size(Rs_SG3, 1), 1)) / (365 * 24 * 60 * 60), 2);

    % Region 4
    Q_RG1(:,5,i) = sum((0.001 .* Rs_RG1(:, MRBidx == 4) .* repmat(area(MRBidx == 4)', size(Rs_RG1, 1), 1)) / (365 * 24 * 60 * 60), 2);
    Q_SG1(:,5,i) = sum((0.001 .* Rs_SG1(:, MRBidx == 4) .* repmat(area(MRBidx == 4)', size(Rs_SG1, 1), 1)) / (365 * 24 * 60 * 60), 2);
    Q_SG2(:,5,i) = sum((0.001 .* Rs_SG2(:, MRBidx == 4) .* repmat(area(MRBidx == 4)', size(Rs_SG2, 1), 1)) / (365 * 24 * 60 * 60), 2);
    Q_SG3(:,5,i) = sum((0.001 .* Rs_SG3(:, MRBidx == 4) .* repmat(area(MRBidx == 4)', size(Rs_SG3, 1), 1)) / (365 * 24 * 60 * 60), 2);

    % Region 5
    Q_RG1(:,6,i) = sum((0.001 .* Rs_RG1(:, MRBidx == 5) .* repmat(area(MRBidx == 5)', size(Rs_RG1, 1), 1)) / (365 * 24 * 60 * 60), 2);
    Q_SG1(:,6,i) = sum((0.001 .* Rs_SG1(:, MRBidx == 5) .* repmat(area(MRBidx == 5)', size(Rs_SG1, 1), 1)) / (365 * 24 * 60 * 60), 2);
    Q_SG2(:,6,i) = sum((0.001 .* Rs_SG2(:, MRBidx == 5) .* repmat(area(MRBidx == 5)', size(Rs_SG2, 1), 1)) / (365 * 24 * 60 * 60), 2);
    Q_SG3(:,6,i) = sum((0.001 .* Rs_SG3(:, MRBidx == 5) .* repmat(area(MRBidx == 5)', size(Rs_SG3, 1), 1)) / (365 * 24 * 60 * 60), 2);

    % Region 6
    Q_RG1(:,7,i) = sum((0.001 .* Rs_RG1(:, MRBidx == 6) .* repmat(area(MRBidx == 6)', size(Rs_RG1, 1), 1)) / (365 * 24 * 60 * 60), 2);
    Q_SG1(:,7,i) = sum((0.001 .* Rs_SG1(:, MRBidx == 6) .* repmat(area(MRBidx == 6)', size(Rs_SG1, 1), 1)) / (365 * 24 * 60 * 60), 2);
    Q_SG2(:,7,i) = sum((0.001 .* Rs_SG2(:, MRBidx == 6) .* repmat(area(MRBidx == 6)', size(Rs_SG2, 1), 1)) / (365 * 24 * 60 * 60), 2);
    Q_SG3(:,7,i) = sum((0.001 .* Rs_SG3(:, MRBidx == 6) .* repmat(area(MRBidx == 6)', size(Rs_SG3, 1), 1)) / (365 * 24 * 60 * 60), 2);
    
end

% time series
subplot(7,2,9:12)
Clim = squeeze(Q_SG1(:,1,:) - Q_RG1(:,1,:));
LULCC = squeeze(Q_SG2(:,1,:) - Q_SG1(:,1,:));
CO2 = squeeze(Q_SG3(:,1,:) - Q_SG2(:,1,:));

% Climate
ci = quantile(Clim, [0.1 0.9], 2);
fill([year(2:end) fliplr(year(2:end))], [ci(2:end,1)' fliplr(ci(2:end,2)')], clr(1,:), 'EdgeColor','none', 'FaceAlpha',0.25)
hold on;
plot(syear:eyear, mean(Clim, 2), '-', 'Color', clr(1,:), 'LineWidth',1)
% LULCC
ci = quantile(LULCC, [0.1 0.9], 2);
fill([year(2:end) fliplr(year(2:end))], [ci(2:end,1)' fliplr(ci(2:end,2)')], clr(2,:), 'EdgeColor','none', 'FaceAlpha',0.25)
plot(syear:eyear, mean(LULCC, 2), '-', 'Color', clr(2,:), 'LineWidth',1)
% CO2
ci = quantile(CO2, [0.1 0.9], 2);
fill([year(2:end) fliplr(year(2:end))], [ci(2:end,1)' fliplr(ci(2:end,2)')], clr(3,:), 'EdgeColor','none', 'FaceAlpha',0.25)
plot(syear:eyear, mean(CO2, 2), '-', 'Color', clr(3,:), 'LineWidth',1)
hold off;
set(gca, 'XLim',[syear eyear], 'YLim',[-4000 6000], 'TickDir','out')
box off;
yl = ylabel('\DeltaQ (m^{3} s^{-1})', 'FontSize',10);
yl.Position(1) = 1892;
hold off;

text(1884,6000,'e','FontSize',12)


% overall
Clim = mean(Q_SG1-Q_RG1, 3);
LULCC = mean(Q_SG2-Q_SG1, 3);
CO2 = mean(Q_SG3-Q_SG2, 3);

dat = [mean(Clim(idx, :))
    mean(LULCC(idx, :))
    mean(CO2(idx, :))];

subplot(7,2,13:14)
b = bar(1:7,dat,'stacked');
b(1).FaceColor = clr(1,:);
b(2).FaceColor = clr(2,:);
b(3).FaceColor = clr(3,:);
ax = gca;
set(ax, 'XLim',[0.5 7.5], 'TickDir','out', 'XTickLabel',{'All','1','2','3','4','5','6'})
box off;
ylb = ylabel(ax, '\DeltaQ (m^{3} s^{-1})', 'FontSize',10);
xlabel('Region', 'FontSize',10)
lgd = legend(ax, '\DeltaQ_{climate}','\DeltaQ_{LULCC}','\DeltaQ_{CO2}', 'Location','north',...
    'Orientation','horizontal');
legend('boxoff')
lgd.Position(2) = 0.44;
ax.Position(2) = 0.07;

text(-0.55, 1100, 'f','FontSize',12)

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/map-preindustrial-contributions-runoff.tif')
close all;

%% plot LULCC vs Climate effects
Clim = squeeze(Q_SG1(:,1,:) - Q_RG1(:,1,:));
LULCC = squeeze(Q_SG2(:,1,:) - Q_SG1(:,1,:));
CO2 = squeeze(Q_SG3(:,1,:) - Q_SG2(:,1,:));

clim_m = mean(Clim, 2);
clim_s = std(Clim, 0, 2) / sqrt(10);

lulcc_m = mean(LULCC, 2);
lulcc_s = std(LULCC, 0, 2) / sqrt(10);

x = [10*round(min(clim_m)/10):10:10*round(max(clim_m)/10)]';
% reduced major axis regression (Trauth, MATLAB Recipes for Earth Sciences)
b1 = std(lulcc_m, "omitnan") / std(clim_m, "omitnan");
b0 = mean(lulcc_m, 'omitnan') - b1*mean(clim_m, 'omitnan');
yhat = x*b1 + b0;
yci = NaN(length(x), 1000);
for i = 1:1000
    ds = datasample([lulcc_m clim_m], length(lulcc_m));
    bs1 = std(ds(:,1), "omitnan") / std(ds(:,2), "omitnan");
    bs0 = mean(ds(:,1), 'omitnan') - b1*mean(ds(:,2), 'omitnan');
    yci(:,i) = x*bs1 + bs0;
end
yci = quantile(yci, [0.025 0.975],2);

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 3.5 3];

fill([x' fliplr(x')], [yci(:,1)' fliplr(yci(:,2)')], [0.8 0.8 0.8], 'EdgeColor','none');
hold on;
plot(x, yhat, 'k-')
scatter(clim_m, lulcc_m, 10, "black", "filled")
set(gca, 'TickDir','out')
box off;
xlabel('\DeltaQ_{Climate} (m^{3} s^{-1})')
ylabel('\DeltaQ_{LULCC} (m^{3} s^{-1})')
xlim = get(gca, 'XLim');
ylim = get(gca, 'YLim');
text(xlim(1)+diff(xlim)*0.02, ylim(2),...
    sprintf('y = %.03fx + %.03f', b1, b0),...
    'FontSize',8)
text(xlim(1)+diff(xlim)*0.02, ylim(2)-0.08*diff(ylim),...
    sprintf('R^{2} = %.02f', corr(lulcc_m, clim_m, 'rows','pairwise')^2),...
    'FontSize',8)

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/climate-vs-lulcc-effect.tif')
close all;

%% Save variables
Q_RG1 = squeeze(Q_RG1(:,1,:));
Q_SG1 = squeeze(Q_SG1(:,1,:));
Q_SG2 = squeeze(Q_SG2(:,1,:));
Q_SG3 = squeeze(Q_SG3(:,1,:));
save('./output/MsTMIP-contributions.mat', 'Clim','LULCC','CO2','Q_SG3','Q_SG2','Q_SG1','Q_RG1','year','models', '-v7.3');


