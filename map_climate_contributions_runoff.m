% map climate contributions runoff

load ./data/McCabeWilliams_WaterBudget;
load ./data/MRB_subregions_GP;
syear = 1981;
eyear = 2010;
[~, ny, nx] = size(Rs_S3_WY);
states = shaperead('usastatehi','UseGeoCoords', true);
latlim=[36 50];
lonlim=[-115 -89];
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
s1 = readmatrix("./modeling/wb.lmrb.awc.runoff.hamonpet.s1a", 'FileType','text');
s2 = readmatrix("./modeling/wb.lmrb.awc.runoff.hamonpet.s2a", 'FileType','text');
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

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 4.75 9];
clr = cbrewer('div','RdBu',20); clr(clr<0) = 0; clr(clr>1) = 1;

subplot(11,4,[9 10 13 14])
axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',4,'MLineLocation',6,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.3,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','MLabelParallel','north',...
        'FLineWidth',1, 'GColor',[0.5 0.5 0.5], 'FontSize',8)
axis off;
axis image;
surfm(lat, lon, squeeze(mean(Rs_S1_WY(year>=syear & year<=eyear, :, :) - Rs_S0_WY(year>=syear & year<=eyear, :, :))));
caxis([-20 20])
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
ttl = title('\DeltaT_{avg}', 'FontSize',12); ttl.Position(2) = 0.84;
ax = gca;
ax.Position(1) = 0.06;
ax.Position(2) = 0.66;
text(-0.18, 0.84, 'c', 'FontSize',12)

subplot(11,4,[1 2 5 6])
axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',4,'MLineLocation',6,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.3,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','MLabelParallel','north',...
        'FLineWidth',1, 'GColor',[0.5 0.5 0.5], 'FontSize',8)
axis off;
axis image;
surfm(lat, lon, squeeze(mean(Rs_S3_WY(year>=syear & year<=eyear, :, :) - Rs_S2_WY(year>=syear & year<=eyear, :, :))));
caxis([-20 20])
colormap(gca, clr)
geoshow(states,'FaceColor','none','EdgeColor',[0.5 0.5 0.5],'LineWidth',0.5)
geoshow(SR1(1,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR2(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR3(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR4(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR5(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR6(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
ttl = title('\DeltaP', 'FontSize',12); ttl.Position(2) = 0.84;
ax = gca;
ax.Position(1) = 0.06;
ax.Position(2) = 0.83;
text(-0.18, 0.84, 'a', 'FontSize',12)

subplot(11,4,[3 4 7 8])
axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',4,'MLineLocation',6,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.3,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','MLabelParallel','north',...
        'FLineWidth',1, 'GColor',[0.5 0.5 0.5], 'FontSize',8)
axis off;
axis image;
surfm(lat, lon, squeeze(mean(Rs_S2_WY(year>=syear & year<=eyear, :, :) - Rs_S1_WY(year>=syear & year<=eyear, :, :))));
caxis([-20 20])
colormap(gca, clr)
geoshow(states,'FaceColor','none','EdgeColor',[0.5 0.5 0.5],'LineWidth',0.5)
geoshow(SR1(1,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR2(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR3(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR4(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR5(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR6(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
ttl = title('\DeltaPET', 'FontSize',12); ttl.Position(2) = 0.84;
ax = gca;
ax.Position(1) = 0.46;
ax.Position(2) = 0.83;
text(-0.18, 0.84, 'b', 'FontSize',12)

subplot(11,4,[11 12 15 16])
axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',4,'MLineLocation',6,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.3,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','MLabelParallel','north',...
        'FLineWidth',1, 'GColor',[0.5 0.5 0.5], 'FontSize',8)
axis off;
axis image;
surfm(lat, lon, squeeze(mean(Rs_S3_WY(year>=syear & year<=eyear, :, :) - Rs_S0_WY(year>=syear & year<=eyear, :, :))));
caxis([-20 20])
colormap(gca, clr)
geoshow(states,'FaceColor','none','EdgeColor',[0.5 0.5 0.5],'LineWidth',0.5)
geoshow(SR1(1,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR2(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR3(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR4(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR5(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
geoshow(SR6(2,1),'FaceColor','none','EdgeColor','k','LineWidth',0.75)
ttl = title('All', 'FontSize',12); ttl.Position(2) = 0.84;
ax = gca;
ax.Position(1) = 0.46;
ax.Position(2) = 0.66;
text(-0.18, 0.84, 'd', 'FontSize',12)

cb = colorbar('eastoutside');
cb.Position = [0.85 0.66 0.03 0.29];
cb.Ticks = -20:2:20;
cb.TickLength = 0.05;
ylabel(cb, '\DeltaR_{s} (mm)', 'FontSize',10)

%% Contributions by season / subbasin
clr = wesanderson('fantasticfox1');
clr = clr(1:3,:);
% clr = [174,1,126;
%     77,146,33;
%     33,102,172]/255;
ysep = 0.085;
ystart = 0.03;

subplot(11,4,41:44)
Q_S0 = reshape(sum((0.001 .* Rs_S0(:, MRBidx == 6) .* repmat(area(MRBidx == 6)', size(Rs_S0, 1), 1)), 2), 12, []) ./ (repmat(days_in_month', 1, length(year)) * 24 * 60 * 60);
Q_S1 = reshape(sum((0.001 .* Rs_S1(:, MRBidx == 6) .* repmat(area(MRBidx == 6)', size(Rs_S1, 1), 1)), 2), 12, []) ./ (repmat(days_in_month', 1, length(year)) * 24 * 60 * 60);
Q_S2 = reshape(sum((0.001 .* Rs_S2(:, MRBidx == 6) .* repmat(area(MRBidx == 6)', size(Rs_S2, 1), 1)), 2), 12, []) ./ (repmat(days_in_month', 1, length(year)) * 24 * 60 * 60);
Q_S3 = reshape(sum((0.001 .* Rs_S3(:, MRBidx == 6) .* repmat(area(MRBidx == 6)', size(Rs_S3, 1), 1)), 2), 12, []) ./ (repmat(days_in_month', 1, length(year)) * 24 * 60 * 60);
dat = [mean(Q_S1(:,year>=syear & year<=eyear), 2)-mean(Q_S0(:,year>=syear & year<=eyear), 2) ...
    mean(Q_S2(:,year>=syear & year<=eyear), 2)-mean(Q_S1(:,year>=syear & year<=eyear), 2) ...
    mean(Q_S3(:,year>=syear & year<=eyear), 2)-mean(Q_S2(:,year>=syear & year<=eyear), 2)];
Q_S0 = filter(bx, a, sum((0.001 .* Rs_S0(:, MRBidx == 6) .* repmat(area(MRBidx == 6)', size(Rs_S0, 1), 1)), 2), [], 1); Q_S0 = Q_S0(mos==9) / (365 * 24 * 60 * 60);
Q_S1 = filter(bx, a, sum((0.001 .* Rs_S1(:, MRBidx == 6) .* repmat(area(MRBidx == 6)', size(Rs_S1, 1), 1)), 2), [], 1); Q_S1 = Q_S1(mos==9) / (365 * 24 * 60 * 60);
Q_S2 = filter(bx, a, sum((0.001 .* Rs_S2(:, MRBidx == 6) .* repmat(area(MRBidx == 6)', size(Rs_S2, 1), 1)), 2), [], 1); Q_S2 = Q_S2(mos==9) / (365 * 24 * 60 * 60);
Q_S3 = filter(bx, a, sum((0.001 .* Rs_S3(:, MRBidx == 6) .* repmat(area(MRBidx == 6)', size(Rs_S3, 1), 1)), 2), [], 1); Q_S3 = Q_S3(mos==9) / (365 * 24 * 60 * 60);
dat(13,:) = [mean(Q_S1(year>=syear & year<=eyear))-mean(Q_S0(year>=syear & year<=eyear)) ...
    mean(Q_S2(year>=syear & year<=eyear))-mean(Q_S1(year>=syear & year<=eyear)) ...
    mean(Q_S3(year>=syear & year<=eyear))-mean(Q_S2(year>=syear & year<=eyear))];
b = bar(1:13,dat,'stacked');
b(1).FaceColor = clr(1,:);
b(2).FaceColor = clr(2,:);
b(3).FaceColor = clr(3,:);
ax = gca;
set(ax, 'XLim',[0.5 13.5], 'TickDir','out', 'XTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D', 'WY'})
box off;
hold on;
scatter(1:13, sum(dat,2), 20, "black","filled")
ax.Position(1) = 0.15;
ax.Position(2) = ystart + 0*ysep;
ylb = ylabel(ax, '\DeltaQ (m^{3} s^{-1})', 'FontSize',8);
ylb.Position(1) = -0.7;
ylim = get(ax, 'YLim');
text(13.5,ylim(2),'Region 6','HorizontalAlignment','right','FontSize',9, 'VerticalAlignment','middle')
text(-1.8,ylim(2),'k','FontSize',12)

subplot(11,4,37:40)
Q_S0 = reshape(sum((0.001 .* Rs_S0(:, MRBidx == 5) .* repmat(area(MRBidx == 5)', size(Rs_S0, 1), 1)), 2), 12, []) ./ (repmat(days_in_month', 1, length(year)) * 24 * 60 * 60);
Q_S1 = reshape(sum((0.001 .* Rs_S1(:, MRBidx == 5) .* repmat(area(MRBidx == 5)', size(Rs_S1, 1), 1)), 2), 12, []) ./ (repmat(days_in_month', 1, length(year)) * 24 * 60 * 60);
Q_S2 = reshape(sum((0.001 .* Rs_S2(:, MRBidx == 5) .* repmat(area(MRBidx == 5)', size(Rs_S2, 1), 1)), 2), 12, []) ./ (repmat(days_in_month', 1, length(year)) * 24 * 60 * 60);
Q_S3 = reshape(sum((0.001 .* Rs_S3(:, MRBidx == 5) .* repmat(area(MRBidx == 5)', size(Rs_S3, 1), 1)), 2), 12, []) ./ (repmat(days_in_month', 1, length(year)) * 24 * 60 * 60);
dat = [mean(Q_S1(:,year>=syear & year<=eyear), 2)-mean(Q_S0(:,year>=syear & year<=eyear), 2) ...
    mean(Q_S2(:,year>=syear & year<=eyear), 2)-mean(Q_S1(:,year>=syear & year<=eyear), 2) ...
    mean(Q_S3(:,year>=syear & year<=eyear), 2)-mean(Q_S2(:,year>=syear & year<=eyear), 2)];
Q_S0 = filter(bx, a, sum((0.001 .* Rs_S0(:, MRBidx == 5) .* repmat(area(MRBidx == 5)', size(Rs_S0, 1), 1)), 2), [], 1); Q_S0 = Q_S0(mos==9) / (365 * 24 * 60 * 60);
Q_S1 = filter(bx, a, sum((0.001 .* Rs_S1(:, MRBidx == 5) .* repmat(area(MRBidx == 5)', size(Rs_S1, 1), 1)), 2), [], 1); Q_S1 = Q_S1(mos==9) / (365 * 24 * 60 * 60);
Q_S2 = filter(bx, a, sum((0.001 .* Rs_S2(:, MRBidx == 5) .* repmat(area(MRBidx == 5)', size(Rs_S2, 1), 1)), 2), [], 1); Q_S2 = Q_S2(mos==9) / (365 * 24 * 60 * 60);
Q_S3 = filter(bx, a, sum((0.001 .* Rs_S3(:, MRBidx == 5) .* repmat(area(MRBidx == 5)', size(Rs_S3, 1), 1)), 2), [], 1); Q_S3 = Q_S3(mos==9) / (365 * 24 * 60 * 60);
dat(13,:) = [mean(Q_S1(year>=syear & year<=eyear))-mean(Q_S0(year>=syear & year<=eyear)) ...
    mean(Q_S2(year>=syear & year<=eyear))-mean(Q_S1(year>=syear & year<=eyear)) ...
    mean(Q_S3(year>=syear & year<=eyear))-mean(Q_S2(year>=syear & year<=eyear))];
b = bar(1:13,dat,'stacked');
b(1).FaceColor = clr(1,:);
b(2).FaceColor = clr(2,:);
b(3).FaceColor = clr(3,:);
ax = gca;
set(ax, 'XLim',[0.5 13.5], 'TickDir','out', 'XTickLabel','')
box off;
hold on;
scatter(1:13, sum(dat,2), 20, "black","filled")
ax.Position(1) = 0.15;
ax.Position(2) = ystart + 1*ysep;
ylb = ylabel(ax, '\DeltaQ (m^{3} s^{-1})', 'FontSize',8);
ylb.Position(1) = -0.7;
ylim = get(ax, 'YLim');
text(13.5,ylim(2),'Region 5','HorizontalAlignment','right','FontSize',9, 'VerticalAlignment','middle')
text(-1.8,ylim(2),'j','FontSize',12)

subplot(11,4,33:36)
Q_S0 = reshape(sum((0.001 .* Rs_S0(:, MRBidx == 4) .* repmat(area(MRBidx == 4)', size(Rs_S0, 1), 1)), 2), 12, []) ./ (repmat(days_in_month', 1, length(year)) * 24 * 60 * 60);
Q_S1 = reshape(sum((0.001 .* Rs_S1(:, MRBidx == 4) .* repmat(area(MRBidx == 4)', size(Rs_S1, 1), 1)), 2), 12, []) ./ (repmat(days_in_month', 1, length(year)) * 24 * 60 * 60);
Q_S2 = reshape(sum((0.001 .* Rs_S2(:, MRBidx == 4) .* repmat(area(MRBidx == 4)', size(Rs_S2, 1), 1)), 2), 12, []) ./ (repmat(days_in_month', 1, length(year)) * 24 * 60 * 60);
Q_S3 = reshape(sum((0.001 .* Rs_S3(:, MRBidx == 4) .* repmat(area(MRBidx == 4)', size(Rs_S3, 1), 1)), 2), 12, []) ./ (repmat(days_in_month', 1, length(year)) * 24 * 60 * 60);
dat = [mean(Q_S1(:,year>=syear & year<=eyear), 2)-mean(Q_S0(:,year>=syear & year<=eyear), 2) ...
    mean(Q_S2(:,year>=syear & year<=eyear), 2)-mean(Q_S1(:,year>=syear & year<=eyear), 2) ...
    mean(Q_S3(:,year>=syear & year<=eyear), 2)-mean(Q_S2(:,year>=syear & year<=eyear), 2)];
Q_S0 = filter(bx, a, sum((0.001 .* Rs_S0(:, MRBidx == 4) .* repmat(area(MRBidx == 4)', size(Rs_S0, 1), 1)), 2), [], 1); Q_S0 = Q_S0(mos==9) / (365 * 24 * 60 * 60);
Q_S1 = filter(bx, a, sum((0.001 .* Rs_S1(:, MRBidx == 4) .* repmat(area(MRBidx == 4)', size(Rs_S1, 1), 1)), 2), [], 1); Q_S1 = Q_S1(mos==9) / (365 * 24 * 60 * 60);
Q_S2 = filter(bx, a, sum((0.001 .* Rs_S2(:, MRBidx == 4) .* repmat(area(MRBidx == 4)', size(Rs_S2, 1), 1)), 2), [], 1); Q_S2 = Q_S2(mos==9) / (365 * 24 * 60 * 60);
Q_S3 = filter(bx, a, sum((0.001 .* Rs_S3(:, MRBidx == 4) .* repmat(area(MRBidx == 4)', size(Rs_S3, 1), 1)), 2), [], 1); Q_S3 = Q_S3(mos==9) / (365 * 24 * 60 * 60);
dat(13,:) = [mean(Q_S1(year>=syear & year<=eyear))-mean(Q_S0(year>=syear & year<=eyear)) ...
    mean(Q_S2(year>=syear & year<=eyear))-mean(Q_S1(year>=syear & year<=eyear)) ...
    mean(Q_S3(year>=syear & year<=eyear))-mean(Q_S2(year>=syear & year<=eyear))];
b = bar(1:13,dat,'stacked');
b(1).FaceColor = clr(1,:);
b(2).FaceColor = clr(2,:);
b(3).FaceColor = clr(3,:);
ax = gca;
set(ax, 'XLim',[0.5 13.5], 'TickDir','out', 'XTickLabel','')
box off;
hold on;
scatter(1:13, sum(dat,2), 20, "black","filled")
ax.Position(1) = 0.15;
ax.Position(2) = ystart + 2*ysep;
ylb = ylabel(ax, '\DeltaQ (m^{3} s^{-1})', 'FontSize',8);
ylb.Position(1) = -0.7;
ylim = get(ax, 'YLim');
text(13.5,ylim(2),'Region 4','HorizontalAlignment','right','FontSize',9, 'VerticalAlignment','middle')
text(-1.8,ylim(2),'i','FontSize',12)

subplot(11,4,29:32)
Q_S0 = reshape(sum((0.001 .* Rs_S0(:, MRBidx == 3) .* repmat(area(MRBidx == 3)', size(Rs_S0, 1), 1)), 2), 12, []) ./ (repmat(days_in_month', 1, length(year)) * 24 * 60 * 60);
Q_S1 = reshape(sum((0.001 .* Rs_S1(:, MRBidx == 3) .* repmat(area(MRBidx == 3)', size(Rs_S1, 1), 1)), 2), 12, []) ./ (repmat(days_in_month', 1, length(year)) * 24 * 60 * 60);
Q_S2 = reshape(sum((0.001 .* Rs_S2(:, MRBidx == 3) .* repmat(area(MRBidx == 3)', size(Rs_S2, 1), 1)), 2), 12, []) ./ (repmat(days_in_month', 1, length(year)) * 24 * 60 * 60);
Q_S3 = reshape(sum((0.001 .* Rs_S3(:, MRBidx == 3) .* repmat(area(MRBidx == 3)', size(Rs_S3, 1), 1)), 2), 12, []) ./ (repmat(days_in_month', 1, length(year)) * 24 * 60 * 60);
dat = [mean(Q_S1(:,year>=syear & year<=eyear), 2)-mean(Q_S0(:,year>=syear & year<=eyear), 2) ...
    mean(Q_S2(:,year>=syear & year<=eyear), 2)-mean(Q_S1(:,year>=syear & year<=eyear), 2) ...
    mean(Q_S3(:,year>=syear & year<=eyear), 2)-mean(Q_S2(:,year>=syear & year<=eyear), 2)];
Q_S0 = filter(bx, a, sum((0.001 .* Rs_S0(:, MRBidx == 3) .* repmat(area(MRBidx == 3)', size(Rs_S0, 1), 1)), 2), [], 1); Q_S0 = Q_S0(mos==9) / (365 * 24 * 60 * 60);
Q_S1 = filter(bx, a, sum((0.001 .* Rs_S1(:, MRBidx == 3) .* repmat(area(MRBidx == 3)', size(Rs_S1, 1), 1)), 2), [], 1); Q_S1 = Q_S1(mos==9) / (365 * 24 * 60 * 60);
Q_S2 = filter(bx, a, sum((0.001 .* Rs_S2(:, MRBidx == 3) .* repmat(area(MRBidx == 3)', size(Rs_S2, 1), 1)), 2), [], 1); Q_S2 = Q_S2(mos==9) / (365 * 24 * 60 * 60);
Q_S3 = filter(bx, a, sum((0.001 .* Rs_S3(:, MRBidx == 3) .* repmat(area(MRBidx == 3)', size(Rs_S3, 1), 1)), 2), [], 1); Q_S3 = Q_S3(mos==9) / (365 * 24 * 60 * 60);
dat(13,:) = [mean(Q_S1(year>=syear & year<=eyear))-mean(Q_S0(year>=syear & year<=eyear)) ...
    mean(Q_S2(year>=syear & year<=eyear))-mean(Q_S1(year>=syear & year<=eyear)) ...
    mean(Q_S3(year>=syear & year<=eyear))-mean(Q_S2(year>=syear & year<=eyear))];
b = bar(1:13,dat,'stacked');
b(1).FaceColor = clr(1,:);
b(2).FaceColor = clr(2,:);
b(3).FaceColor = clr(3,:);
ax = gca;
set(ax, 'XLim',[0.5 13.5], 'TickDir','out', 'XTickLabel','')
box off;
hold on;
scatter(1:13, sum(dat,2), 20, "black","filled")
ax.Position(1) = 0.15;
ax.Position(2) = ystart + 3*ysep;
ylb = ylabel(ax, '\DeltaQ (m^{3} s^{-1})', 'FontSize',8);
ylb.Position(1) = -0.7;
ylim = get(ax, 'YLim');
text(13.5,ylim(2),'Region 3','HorizontalAlignment','right','FontSize',9, 'VerticalAlignment','middle')
text(-1.8,ylim(2),'h','FontSize',12)

subplot(11,4,25:28)
Q_S0 = reshape(sum((0.001 .* Rs_S0(:, MRBidx == 2) .* repmat(area(MRBidx == 2)', size(Rs_S0, 1), 1)), 2), 12, []) ./ (repmat(days_in_month', 1, length(year)) * 24 * 60 * 60);
Q_S1 = reshape(sum((0.001 .* Rs_S1(:, MRBidx == 2) .* repmat(area(MRBidx == 2)', size(Rs_S1, 1), 1)), 2), 12, []) ./ (repmat(days_in_month', 1, length(year)) * 24 * 60 * 60);
Q_S2 = reshape(sum((0.001 .* Rs_S2(:, MRBidx == 2) .* repmat(area(MRBidx == 2)', size(Rs_S2, 1), 1)), 2), 12, []) ./ (repmat(days_in_month', 1, length(year)) * 24 * 60 * 60);
Q_S3 = reshape(sum((0.001 .* Rs_S3(:, MRBidx == 2) .* repmat(area(MRBidx == 2)', size(Rs_S3, 1), 1)), 2), 12, []) ./ (repmat(days_in_month', 1, length(year)) * 24 * 60 * 60);
dat = [mean(Q_S1(:,year>=syear & year<=eyear), 2)-mean(Q_S0(:,year>=syear & year<=eyear), 2) ...
    mean(Q_S2(:,year>=syear & year<=eyear), 2)-mean(Q_S1(:,year>=syear & year<=eyear), 2) ...
    mean(Q_S3(:,year>=syear & year<=eyear), 2)-mean(Q_S2(:,year>=syear & year<=eyear), 2)];
Q_S0 = filter(bx, a, sum((0.001 .* Rs_S0(:, MRBidx == 2) .* repmat(area(MRBidx == 2)', size(Rs_S0, 1), 1)), 2), [], 1); Q_S0 = Q_S0(mos==9) / (365 * 24 * 60 * 60);
Q_S1 = filter(bx, a, sum((0.001 .* Rs_S1(:, MRBidx == 2) .* repmat(area(MRBidx == 2)', size(Rs_S1, 1), 1)), 2), [], 1); Q_S1 = Q_S1(mos==9) / (365 * 24 * 60 * 60);
Q_S2 = filter(bx, a, sum((0.001 .* Rs_S2(:, MRBidx == 2) .* repmat(area(MRBidx == 2)', size(Rs_S2, 1), 1)), 2), [], 1); Q_S2 = Q_S2(mos==9) / (365 * 24 * 60 * 60);
Q_S3 = filter(bx, a, sum((0.001 .* Rs_S3(:, MRBidx == 2) .* repmat(area(MRBidx == 2)', size(Rs_S3, 1), 1)), 2), [], 1); Q_S3 = Q_S3(mos==9) / (365 * 24 * 60 * 60);
dat(13,:) = [mean(Q_S1(year>=syear & year<=eyear))-mean(Q_S0(year>=syear & year<=eyear)) ...
    mean(Q_S2(year>=syear & year<=eyear))-mean(Q_S1(year>=syear & year<=eyear)) ...
    mean(Q_S3(year>=syear & year<=eyear))-mean(Q_S2(year>=syear & year<=eyear))];
b = bar(1:13,dat,'stacked');
b(1).FaceColor = clr(1,:);
b(2).FaceColor = clr(2,:);
b(3).FaceColor = clr(3,:);
ax = gca;
set(ax, 'XLim',[0.5 13.5], 'TickDir','out', 'XTickLabel','')
box off;
hold on;
scatter(1:13, sum(dat,2), 20, "black","filled")
ax.Position(1) = 0.15;
ax.Position(2) = ystart + 4*ysep;
ylb = ylabel(ax, '\DeltaQ (m^{3} s^{-1})', 'FontSize',8);
ylb.Position(1) = -0.7;
ylim = get(ax, 'YLim');
text(13.5,ylim(2),'Region 2','HorizontalAlignment','right','FontSize',9, 'VerticalAlignment','middle')
text(-1.8,ylim(2),'g','FontSize',12)

subplot(11,4,21:24)
Q_S0 = reshape(sum((0.001 .* Rs_S0(:, MRBidx == 1) .* repmat(area(MRBidx == 1)', size(Rs_S0, 1), 1)), 2), 12, []) ./ (repmat(days_in_month', 1, length(year)) * 24 * 60 * 60);
Q_S1 = reshape(sum((0.001 .* Rs_S1(:, MRBidx == 1) .* repmat(area(MRBidx == 1)', size(Rs_S1, 1), 1)), 2), 12, []) ./ (repmat(days_in_month', 1, length(year)) * 24 * 60 * 60);
Q_S2 = reshape(sum((0.001 .* Rs_S2(:, MRBidx == 1) .* repmat(area(MRBidx == 1)', size(Rs_S2, 1), 1)), 2), 12, []) ./ (repmat(days_in_month', 1, length(year)) * 24 * 60 * 60);
Q_S3 = reshape(sum((0.001 .* Rs_S3(:, MRBidx == 1) .* repmat(area(MRBidx == 1)', size(Rs_S3, 1), 1)), 2), 12, []) ./ (repmat(days_in_month', 1, length(year)) * 24 * 60 * 60);
dat = [mean(Q_S1(:,year>=syear & year<=eyear), 2)-mean(Q_S0(:,year>=syear & year<=eyear), 2) ...
    mean(Q_S2(:,year>=syear & year<=eyear), 2)-mean(Q_S1(:,year>=syear & year<=eyear), 2) ...
    mean(Q_S3(:,year>=syear & year<=eyear), 2)-mean(Q_S2(:,year>=syear & year<=eyear), 2)];
Q_S0 = filter(bx, a, sum((0.001 .* Rs_S0(:, MRBidx == 1) .* repmat(area(MRBidx == 1)', size(Rs_S0, 1), 1)), 2), [], 1); Q_S0 = Q_S0(mos==9) / (365 * 24 * 60 * 60);
Q_S1 = filter(bx, a, sum((0.001 .* Rs_S1(:, MRBidx == 1) .* repmat(area(MRBidx == 1)', size(Rs_S1, 1), 1)), 2), [], 1); Q_S1 = Q_S1(mos==9) / (365 * 24 * 60 * 60);
Q_S2 = filter(bx, a, sum((0.001 .* Rs_S2(:, MRBidx == 1) .* repmat(area(MRBidx == 1)', size(Rs_S2, 1), 1)), 2), [], 1); Q_S2 = Q_S2(mos==9) / (365 * 24 * 60 * 60);
Q_S3 = filter(bx, a, sum((0.001 .* Rs_S3(:, MRBidx == 1) .* repmat(area(MRBidx == 1)', size(Rs_S3, 1), 1)), 2), [], 1); Q_S3 = Q_S3(mos==9) / (365 * 24 * 60 * 60);
dat(13,:) = [mean(Q_S1(year>=syear & year<=eyear))-mean(Q_S0(year>=syear & year<=eyear)) ...
    mean(Q_S2(year>=syear & year<=eyear))-mean(Q_S1(year>=syear & year<=eyear)) ...
    mean(Q_S3(year>=syear & year<=eyear))-mean(Q_S2(year>=syear & year<=eyear))];
b = bar(1:13,dat,'stacked');
b(1).FaceColor = clr(1,:);
b(2).FaceColor = clr(2,:);
b(3).FaceColor = clr(3,:);
ax = gca;
set(ax, 'XLim',[0.5 13.5], 'TickDir','out', 'XTickLabel','')
box off;
hold on;
scatter(1:13, sum(dat,2), 20, "black","filled")
ax.Position(1) = 0.15;
ax.Position(2) = ystart + 5*ysep;
ylb = ylabel(ax, '\DeltaQ (m^{3} s^{-1})', 'FontSize',8);
ylb.Position(1) = -0.7;
ylim = get(ax, 'YLim');
text(13.5,ylim(2),'Region 1','HorizontalAlignment','right','FontSize',9, 'VerticalAlignment','middle')
text(-1.8,ylim(2),'f','FontSize',12)

subplot(11,4,17:20)
Q_S0 = reshape(sum((0.001 .* Rs_S0(:, MRBidx > 0) .* repmat(area(MRBidx > 0)', size(Rs_S0, 1), 1)), 2), 12, []) ./ (repmat(days_in_month', 1, length(year)) * 24 * 60 * 60);
Q_S1 = reshape(sum((0.001 .* Rs_S1(:, MRBidx > 0) .* repmat(area(MRBidx > 0)', size(Rs_S1, 1), 1)), 2), 12, []) ./ (repmat(days_in_month', 1, length(year)) * 24 * 60 * 60);
Q_S2 = reshape(sum((0.001 .* Rs_S2(:, MRBidx > 0) .* repmat(area(MRBidx > 0)', size(Rs_S2, 1), 1)), 2), 12, []) ./ (repmat(days_in_month', 1, length(year)) * 24 * 60 * 60);
Q_S3 = reshape(sum((0.001 .* Rs_S3(:, MRBidx > 0) .* repmat(area(MRBidx > 0)', size(Rs_S3, 1), 1)), 2), 12, []) ./ (repmat(days_in_month', 1, length(year)) * 24 * 60 * 60);
dat = [mean(Q_S1(:,year>=syear & year<=eyear), 2)-mean(Q_S0(:,year>=syear & year<=eyear), 2) ...
    mean(Q_S2(:,year>=syear & year<=eyear), 2)-mean(Q_S1(:,year>=syear & year<=eyear), 2) ...
    mean(Q_S3(:,year>=syear & year<=eyear), 2)-mean(Q_S2(:,year>=syear & year<=eyear), 2)];
Q_S0 = filter(bx, a, sum((0.001 .* Rs_S0(:, MRBidx > 0) .* repmat(area(MRBidx > 0)', size(Rs_S0, 1), 1)), 2), [], 1); Q_S0 = Q_S0(mos==9) / (365 * 24 * 60 * 60);
Q_S1 = filter(bx, a, sum((0.001 .* Rs_S1(:, MRBidx > 0) .* repmat(area(MRBidx > 0)', size(Rs_S1, 1), 1)), 2), [], 1); Q_S1 = Q_S1(mos==9) / (365 * 24 * 60 * 60);
Q_S2 = filter(bx, a, sum((0.001 .* Rs_S2(:, MRBidx > 0) .* repmat(area(MRBidx > 0)', size(Rs_S2, 1), 1)), 2), [], 1); Q_S2 = Q_S2(mos==9) / (365 * 24 * 60 * 60);
Q_S3 = filter(bx, a, sum((0.001 .* Rs_S3(:, MRBidx > 0) .* repmat(area(MRBidx > 0)', size(Rs_S3, 1), 1)), 2), [], 1); Q_S3 = Q_S3(mos==9) / (365 * 24 * 60 * 60);
dat(13,:) = [mean(Q_S1(year>=syear & year<=eyear))-mean(Q_S0(year>=syear & year<=eyear)) ...
    mean(Q_S2(year>=syear & year<=eyear))-mean(Q_S1(year>=syear & year<=eyear)) ...
    mean(Q_S3(year>=syear & year<=eyear))-mean(Q_S2(year>=syear & year<=eyear))];
b = bar(1:13,dat,'stacked');
b(1).FaceColor = clr(1,:);
b(2).FaceColor = clr(2,:);
b(3).FaceColor = clr(3,:);
ax = gca;
set(ax, 'XLim',[0.5 13.5], 'TickDir','out', 'XTickLabel','')
box off;
hold on;
scatter(1:13, sum(dat,2), 20, "black","filled")
ax.Position(1) = 0.15;
ax.Position(2) = ystart + 6*ysep;
ylb = ylabel(ax, '\DeltaQ (m^{3} s^{-1})', 'FontSize',8);
ylb.Position(1) = -0.7;
ylim = get(ax, 'YLim');
text(13.5,ylim(2),'All regions','HorizontalAlignment','right','FontSize',9, 'VerticalAlignment','middle')
lgd = legend('\DeltaT_{avg}','\DeltaPET','\DeltaP', 'Location','northoutside', 'Orientation', 'horizontal');
lgd.Position(2) = 0.61;
legend('boxoff')
text(-1.8,ylim(2),'e','FontSize',12)

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300',['./output/map-climate-contributions-runoff-',num2str(syear),'-',num2str(eyear),'.tif'])
close all;


