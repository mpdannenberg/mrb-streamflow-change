load ./data/McCabeWilliams_WaterBudget;
load ./data/MRB_subregions_GP;
syear = 1981;
eyear = 2010;
[~, ny, nx] = size(Rs_S3_WY);
% states = shaperead('usastatehi','UseGeoCoords', true);
% latlim=[36 50];
% lonlim=[-115 -89];
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


%% Contributions by season / subbasin
clr = wesanderson('fantasticfox1');
clr = clr(1:3,:);

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 0 4.75 6];
ysep = 0.16;
ystart = 0.05;

subplot(6,4,21:24)
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
text(-1.8,ylim(2),'f','FontSize',12)

subplot(6,4,17:20)
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
text(-1.8,ylim(2),'e','FontSize',12)

subplot(6,4,13:16)
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
text(-1.8,ylim(2),'d','FontSize',12)

subplot(6,4,9:12)
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
text(-1.8,ylim(2),'c','FontSize',12)

subplot(6,4,5:8)
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
text(-1.8,ylim(2),'b','FontSize',12)

subplot(6,4,1:4)
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
text(-1.8,ylim(2),'a','FontSize',12)

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300',['./output/monthly-climate-contributions-runoff-',num2str(syear),'-',num2str(eyear),'.tif'])
close all;

