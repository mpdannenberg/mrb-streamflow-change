% Calculate runoff ratio and trend for each dataset

load ./data/MRB_subregions_GP;
states = shaperead('usastatehi','UseGeoCoords', true);
latlim=[36 50];
lonlim=[-115 -89];
windowSize = 12;
bx = ones(1,windowSize);
a = 1;

%% Water budget model
% get subbasin number
load ./data/McCabeWilliams_WaterBudget;
[~, ny, nx] = size(Rs_S3_WY);
cellsize = 0.25;

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

e = referenceEllipsoid('World Geodetic System 1984');
area = areaquad(reshape(LAT-(cellsize/2),[],1),reshape(LON-(cellsize/2),[],1),reshape(LAT+(cellsize/2),[],1),reshape(LON+(cellsize/2),[],1),e);
area = reshape(area, ny, nx); 

clear LAT LON LatLon IN ON idx e;

% Organize precip data and calculate runoff ratio
clat = readmatrix('./modeling/lat.txt');
clon = readmatrix('./modeling/lon.txt');
clim = readmatrix("./modeling/prec.csv");
yr = clim(:,1);
mo = clim(:,2);
P = NaN(size(clim,1), ny, nx);
for i = 1:length(clat)
    xi = find(lon == clon(i));
    yi = find(lat == clat(i));
    P(:, yi, xi) = clim(:,i+2);
end
P12 = filter(bx, a, P, [], 1);
P12(1:(windowSize-1), :, :) = NaN;
WYP = P12(mo==9, :, :);
Ptot = sum((0.001 * WYP(:, MRBidx>0)) .* repmat(area(MRBidx>0)', length(year),1),2) / (365*24*60*60);
RR_WBM = Q_S3_WY ./ Ptot;
Ptot_WBM = Ptot; year_WBM = year;
clear a bx clat clon i mo P P12 WYP windowSize xi yi yr Ptot year;

%% gage runoff ratio
T = readtable('./data/Missouri_Big6_UpdatedUnregulatedFlow_DailyMonthlyWY_avgCMS.xlsx',...
    'Sheet','WaterYear_cms');
T.Properties.VariableNames = {'Year','FTRA','SUX','NCNE','MKC','MNMO','HEMO'};
[~, ia, ib] = intersect(year_WBM, T.Year);
RR_gage = T.HEMO(ib) ./ Ptot_WBM(ia);
year_gage = T.Year(ib);

%% MsTMIP runoff ratio
load ./data/MsTMIP_WaterBudget;
[nt,ny,nx,nk] = size(ET_SG3);
cellsize = 0.5;

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

e = referenceEllipsoid('World Geodetic System 1984');
area = areaquad(reshape(LAT-(cellsize/2),[],1),reshape(LON-(cellsize/2),[],1),reshape(LAT+(cellsize/2),[],1),reshape(LON+(cellsize/2),[],1),e);
area = reshape(area, ny, nx); 

clear LAT LON LatLon IN ON idx e;

temp = NaN(1,length(area(MRBidx>0))); temp(:,:) = area(MRBidx>0);
Ptot = sum(P(:,MRBidx>0) .* repmat(temp, nt, 1), 2) * (1/365) * (1/24) * (1/60) * (1/60);
WB_SG3 = permute(repmat(P, 1,1,1,nk) - ET_SG3, [1 4 2 3]);
temp = NaN(1,1,length(area(MRBidx>0))); temp(:,:,:) = area(MRBidx>0);
Q = sum(WB_SG3(:,:,MRBidx>0) .* repmat(temp, nt, nk, 1), 3) * (1/365) * (1/24) * (1/60) * (1/60);
RR_mstmip = Q ./ repmat(Ptot,1,nk);
Ptot_mstmip = Ptot; year_mstmip = year; clear Ptot year;

%% Make plot
clr = wesanderson('fantasticfox1');
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6.5 3];

ci = quantile(RR_mstmip, [0.1 0.9],2);

fill([year_mstmip(2:end) fliplr(year_mstmip(2:end))],...
    [ci(2:end,1)' fliplr(ci(2:end,2)')],...
    clr(1,:), 'EdgeColor','none', 'FaceAlpha',0.3)
hold on;
p2 = plot(year_mstmip', mean(RR_mstmip, 2), '-', 'Color',clr(1,:), 'LineWidth',1.2);
p3 = plot(year_WBM', RR_WBM, '-', 'Color',clr(2,:), 'LineWidth',1.2);
p1 = plot(year_gage, RR_gage, 'k-', 'LineWidth',1.2);
set(gca, 'YLim',[0 0.6], 'XLim',[1930 2010], 'TickDir','out')
box off;
ylabel('Runoff ratio (Q/P)')
ylim = get(gca, 'YLim');
xlabel('Water year')

mdl = fitlm(year_gage(year_gage >= 1931 & year_gage <= 2010),...
    RR_gage(year_gage >= 1931 & year_gage <= 2010));
text(1932, ylim(2), ['{\bfGage}: ', sprintf('%0.03f ',mdl.Coefficients.Estimate(2)*10), 'decade^{-1} (p = ',sprintf('%0.02f',mdl.Coefficients.pValue(2)),')'], 'FontSize',8)
mdl = fitlm(year_mstmip(year_mstmip >= 1931 & year_mstmip <= 2010),...
    mean(RR_mstmip(year_mstmip >= 1931 & year_mstmip <= 2010,:),2));
text(1932, ylim(2)-0.07*diff(ylim), ['{\bfMsTMIP}: ', sprintf('%0.03f ',mdl.Coefficients.Estimate(2)*10), 'decade^{-1} (p = ',sprintf('%0.02f',mdl.Coefficients.pValue(2)),')'], 'FontSize',8, 'Color', clr(1,:))
mdl = fitlm(year_WBM(year_WBM >= 1931 & year_WBM <= 2010),...
    RR_WBM(year_WBM >= 1931 & year_WBM <= 2010));
text(1932, ylim(2)-0.14*diff(ylim), ['{\bfMW11}: ', sprintf('%0.03f ',mdl.Coefficients.Estimate(2)*10), 'decade^{-1} (p = ',sprintf('%0.02f',mdl.Coefficients.pValue(2)),')'], 'FontSize',8, 'Color',clr(2,:))

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/mrb-runoff-ratios.tif')
close all;

