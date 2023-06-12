
T = readtable('./data/Missouri_Big6_UpdatedUnregulatedFlow_DailyMonthlyWY_avgCMS.xlsx',...
    'Sheet','WaterYear_cms');
T.Properties.VariableNames = {'Year','FTRA','SUX','NCNE','MKC','MNMO','HEMO'};

basins = shaperead('./data/LMRBgagebasins/lmrb_shapefiles.shp');

%% Find MsTMIP grid cells for each basin and calculate area
load ./data/MsTMIP_WaterBudget;
[nt,ny,nx,nk] = size(ET_SG3);

[LON, LAT] = meshgrid(lon, lat);
LatLon = [reshape(LAT, [], 1) reshape(LON, [], 1)];

[IN, ON] = inpolygon(LatLon(:,1), LatLon(:,2), basins(1).Y, basins(1).X);
idx = IN | ON;
FTRA = reshape(idx, ny, nx);

[IN, ON] = inpolygon(LatLon(:,1), LatLon(:,2), basins(2).Y, basins(2).X);
idx = IN | ON;
SUX = reshape(idx, ny, nx);

[IN, ON] = inpolygon(LatLon(:,1), LatLon(:,2), basins(3).Y, basins(3).X);
idx = IN | ON;
NCNE = reshape(idx, ny, nx);

[IN, ON] = inpolygon(LatLon(:,1), LatLon(:,2), basins(4).Y, basins(4).X);
idx = IN | ON;
MKC = reshape(idx, ny, nx);

[IN, ON] = inpolygon(LatLon(:,1), LatLon(:,2), basins(5).Y, basins(5).X);
idx = IN | ON;
MNMO = reshape(idx, ny, nx);

[IN, ON] = inpolygon(LatLon(:,1), LatLon(:,2), basins(6).Y, basins(6).X);
idx = IN | ON;
HEMO = reshape(idx, ny, nx);

e = referenceEllipsoid('World Geodetic System 1984');
area = areaquad(reshape(LAT-0.25,[],1),reshape(LON-0.25,[],1),reshape(LAT+0.25,[],1),reshape(LON+0.25,[],1),e);
area = reshape(area, ny, nx); 

clear LAT LON LatLon IN ON idx e;

%% Calculate total streamflow in each basin
WB_RG1 = permute(repmat(mean(P(1:30,:,:),'omitnan'),nt,1,1,nk) - ET_RG1, [1 4 2 3]);
WB_SG1 = permute(repmat(P, 1,1,1,nk) - ET_SG1, [1 4 2 3]);
WB_SG2 = permute(repmat(P, 1,1,1,nk) - ET_SG2, [1 4 2 3]);
WB_SG3 = permute(repmat(P, 1,1,1,nk) - ET_SG3, [1 4 2 3]);

%% Make figure showing comparison between gage data and MsTMIP
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6.5 6];

[~, ia, ib] = intersect(year, T.Year);

% FTRA
temp = NaN(1,1,length(area(FTRA))); temp(:,:,:) = area(FTRA);
Q = sum(WB_SG3(:,:,FTRA) .* repmat(temp, nt, nk, 1), 3) * (1/365) * (1/24) * (1/60) * (1/60);

subplot(6,4,[1 3])
plot(year(ia)', Q(ia,:), '-', 'Color',[0.8 0.8 0.8])
hold on;
p1 = plot(year(ia)', mean(Q(ia,:), 2), '-', 'Color',[0.5 0.5 0.5], 'LineWidth',1.5);
p2 = plot(T.Year(ib), T.FTRA(ib), 'k-', 'LineWidth',1.5);
set(gca, 'TickDir','out', 'XTickLabel','')
box off;
ylabel('Q (m^{3} s^{-1})')
ylim = get(gca, 'YLim');
text(1931, ylim(2), 'Fort Randall Dam', 'FontSize',10, 'VerticalAlignment','middle')
lgd = legend([p2 p1], 'Gage', 'MsTMIP', 'Location','north', 'Orientation','horizontal');
lgd.Position(1) = 0.45;
lgd.Position(2) = 0.92;
legend('boxoff')

subplot(6,4,4)
lim = [0 3200];
%plot(lim, lim, '-', 'Color',[0.5 0.5 0.5])
hold on;
scatter(T.FTRA(ib), mean(Q(ia,:), 2), 10, "black","filled")
%set(gca, 'TickDir','out', 'YLim',lim, 'XLim',lim, 'TickLength',[0.04 0])
box off;
ax = gca;
ax.Position(1) = 0.82;
ax.Position(3) = 0.12;
ylabel('MsTMIP')
r = corr(T.FTRA(ib), mean(Q(ia,:), 2), 'rows','pairwise');
text(lim(2), lim(1), sprintf('R = %.02f', r), 'HorizontalAlignment','right', 'VerticalAlignment','bottom', 'FontSize',7)

% SUX
temp = NaN(1,1,length(area(SUX))); temp(:,:,:) = area(SUX);
Q = sum(WB_SG3(:,:,SUX) .* repmat(temp, nt, nk, 1), 3) * (1/365) * (1/24) * (1/60) * (1/60);

subplot(6,4,[5 7])
plot(year(ia)', Q(ia,:), '-', 'Color',[0.8 0.8 0.8])
hold on;
plot(year(ia)', mean(Q(ia,:), 2), '-', 'Color',[0.5 0.5 0.5], 'LineWidth',1.5)
plot(T.Year(ib), T.SUX(ib), 'k-', 'LineWidth',1.5)
set(gca, 'TickDir','out', 'XTickLabel','')
box off;
ylabel('Q (m^{3} s^{-1})')
ylim = get(gca, 'YLim');
text(1931, ylim(2), 'Sioux City, IA', 'FontSize',10, 'VerticalAlignment','middle')

subplot(6,4,8)
lim = [0 4000];
%plot(lim, lim, '-', 'Color',[0.5 0.5 0.5])
hold on;
scatter(T.SUX(ib), mean(Q(ia,:), 2), 10, "black","filled")
%set(gca, 'TickDir','out', 'YLim',lim, 'XLim',lim, 'TickLength',[0.04 0])
box off;
ax = gca;
ax.Position(1) = 0.82;
ax.Position(3) = 0.12;
ylabel('MsTMIP')
r = corr(T.SUX(ib), mean(Q(ia,:), 2), 'rows','pairwise');
text(lim(2), lim(1), sprintf('R = %.02f', r), 'HorizontalAlignment','right', 'VerticalAlignment','bottom', 'FontSize',7)

% NCNE
temp = NaN(1,1,length(area(NCNE))); temp(:,:,:) = area(NCNE);
Q = sum(WB_SG3(:,:,NCNE) .* repmat(temp, nt, nk, 1), 3) * (1/365) * (1/24) * (1/60) * (1/60);

subplot(6,4,[9 11])
plot(year(ia)', Q(ia,:), '-', 'Color',[0.8 0.8 0.8])
hold on;
plot(year(ia)', mean(Q(ia,:), 2), '-', 'Color',[0.5 0.5 0.5], 'LineWidth',1.5)
plot(T.Year(ib), T.NCNE(ib), 'k-', 'LineWidth',1.5)
set(gca, 'TickDir','out', 'XTickLabel','')
box off;
ylabel('Q (m^{3} s^{-1})')
ylim = get(gca, 'YLim');
text(1931, ylim(2), 'Nebraska City, NE', 'FontSize',10, 'VerticalAlignment','middle')

subplot(6,4,12)
lim = [0 5400];
%plot(lim, lim, '-', 'Color',[0.5 0.5 0.5])
hold on;
scatter(T.NCNE(ib), mean(Q(ia,:), 2), 10, "black","filled")
%set(gca, 'TickDir','out', 'YLim',lim, 'XLim',lim, 'TickLength',[0.04 0])
box off;
ax = gca;
ax.Position(1) = 0.82;
ax.Position(3) = 0.12;
ylabel('MsTMIP')
r = corr(T.NCNE(ib), mean(Q(ia,:), 2), 'rows','pairwise');
text(lim(2), lim(1), sprintf('R = %.02f', r), 'HorizontalAlignment','right', 'VerticalAlignment','bottom', 'FontSize',7)

% MKC
temp = NaN(1,1,length(area(MKC))); temp(:,:,:) = area(MKC);
Q = sum(WB_SG3(:,:,MKC) .* repmat(temp, nt, nk, 1), 3) * (1/365) * (1/24) * (1/60) * (1/60);

subplot(6,4,[13 15])
plot(year(ia)', Q(ia,:), '-', 'Color',[0.8 0.8 0.8])
hold on;
plot(year(ia)', mean(Q(ia,:), 2), '-', 'Color',[0.5 0.5 0.5], 'LineWidth',1.5)
plot(T.Year(ib), T.MKC(ib), 'k-', 'LineWidth',1.5)
set(gca, 'TickDir','out', 'XTickLabel','')
box off;
ylabel('Q (m^{3} s^{-1})')
ylim = get(gca, 'YLim');
text(1931, ylim(2), 'Kansas City, MO', 'FontSize',10, 'VerticalAlignment','middle')

subplot(6,4,16)
lim = [0 8000];
%plot(lim, lim, '-', 'Color',[0.5 0.5 0.5])
hold on;
scatter(T.MKC(ib), mean(Q(ia,:), 2), 10, "black","filled")
%set(gca, 'TickDir','out', 'YLim',lim, 'XLim',lim, 'TickLength',[0.04 0])
box off;
ax = gca;
ax.Position(1) = 0.82;
ax.Position(3) = 0.12;
ylabel('MsTMIP')
r = corr(T.MKC(ib), mean(Q(ia,:), 2), 'rows','pairwise');
text(lim(2), lim(1), sprintf('R = %.02f', r), 'HorizontalAlignment','right', 'VerticalAlignment','bottom', 'FontSize',7)

% MNMO
temp = NaN(1,1,length(area(MNMO))); temp(:,:,:) = area(MNMO);
Q = sum(WB_SG3(:,:,MNMO) .* repmat(temp, nt, nk, 1), 3) * (1/365) * (1/24) * (1/60) * (1/60);

subplot(6,4,[17 19])
plot(year(ia)', Q(ia,:), '-', 'Color',[0.8 0.8 0.8])
hold on;
plot(year(ia)', mean(Q(ia,:), 2), '-', 'Color',[0.5 0.5 0.5], 'LineWidth',1.5)
plot(T.Year(ib), T.MNMO(ib), 'k-', 'LineWidth',1.5)
set(gca, 'TickDir','out', 'XTickLabel','')
box off;
ylabel('Q (m^{3} s^{-1})')
ylim = get(gca, 'YLim');
text(1931, ylim(2), 'Boonville, MO', 'FontSize',10, 'VerticalAlignment','middle')

subplot(6,4,20)
lim = [0 9000];
%plot(lim, lim, '-', 'Color',[0.5 0.5 0.5])
hold on;
scatter(T.MNMO(ib), mean(Q(ia,:), 2), 10, "black","filled")
%set(gca, 'TickDir','out', 'YLim',lim, 'XLim',lim, 'TickLength',[0.04 0])
box off;
ax = gca;
ax.Position(1) = 0.82;
ax.Position(3) = 0.12;
ylabel('MsTMIP')
r = corr(T.MNMO(ib), mean(Q(ia,:), 2), 'rows','pairwise');
text(lim(2), lim(1), sprintf('R = %.02f', r), 'HorizontalAlignment','right', 'VerticalAlignment','bottom', 'FontSize',7)

% HEMO
temp = NaN(1,1,length(area(HEMO))); temp(:,:,:) = area(HEMO);
Q = sum(WB_SG3(:,:,HEMO) .* repmat(temp, nt, nk, 1), 3) * (1/365) * (1/24) * (1/60) * (1/60);

subplot(6,4,[21 23])
plot(year(ia)', Q(ia,:), '-', 'Color',[0.8 0.8 0.8])
hold on;
plot(year(ia)', mean(Q(ia,:), 2), '-', 'Color',[0.5 0.5 0.5], 'LineWidth',1.5)
plot(T.Year(ib), T.HEMO(ib), 'k-', 'LineWidth',1.5)
set(gca, 'TickDir','out')
box off;
ylabel('Q (m^{3} s^{-1})')
ylim = get(gca, 'YLim');
text(1931, ylim(2), 'Hermann, MO', 'FontSize',10, 'VerticalAlignment','middle')

subplot(6,4,24)
lim = [0 10500];
%plot(lim, lim, '-', 'Color',[0.5 0.5 0.5])
hold on;
scatter(T.HEMO(ib), mean(Q(ia,:), 2), 10, "black","filled")
%set(gca, 'TickDir','out', 'YLim',lim, 'XLim',lim, 'TickLength',[0.04 0])
box off;
ax = gca;
ax.Position(1) = 0.82;
ax.Position(3) = 0.12;
ylabel('MsTMIP')
xlabel('Gage')
r = corr(T.HEMO(ib), mean(Q(ia,:), 2), 'rows','pairwise');
text(lim(2), lim(1), sprintf('R = %.02f', r), 'HorizontalAlignment','right', 'VerticalAlignment','bottom', 'FontSize',7)

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r600','./output/MsTMIP-vs-gages.tif')
close all;


