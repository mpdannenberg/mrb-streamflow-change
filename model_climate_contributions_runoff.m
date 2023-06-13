% Test tree-ring reconstruction of WY P-E in LMR
cd modeling;
cellsize = 0.25; % grid resolution in degrees
year = 1901:2021;

%% Calculate water-year climatic water balance (CWB) and aridity index (AI)
clatlim = [32 50]; lat = (clatlim(1) + 0.25/2):0.25:(clatlim(2) - 0.25/2);
clonlim = [-115 -85]; lon = (clonlim(1) + 0.25/2):0.25:(clonlim(2) - 0.25/2);
yr = reshape(repmat(1901:2021, 12, 1), [], 1);
mos = repmat(1:12, 1, length(1901:2021))';
clat = readmatrix('lat.txt');
clon = readmatrix('lon.txt');
nt = length(yr); ny = length(lat); nx = length(lon);

s0 = readmatrix("wb.lmrb.awc.runoff.hamonpet.s0", 'FileType','text');
s1 = readmatrix("wb.lmrb.awc.runoff.hamonpet.s1", 'FileType','text');
s2 = readmatrix("wb.lmrb.awc.runoff.hamonpet.s2", 'FileType','text');
s3 = readmatrix("wb.lmrb.awc.runoff.hamonpet.s3", 'FileType','text');

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

% Find grid cells in MRB and calculate area of each cell
load ../data/MRB_subregions_GP;
[LON, LAT] = meshgrid(lon, lat);
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
area = areaquad(reshape(LAT-(cellsize/2),[],1),reshape(LON-(cellsize/2),[],1),reshape(LAT+(cellsize/2),[],1),reshape(LON+(cellsize/2),[],1),e);
area = reshape(area, ny, nx); 

clear SR* LAT LON LatLon IN ON idx;

% calculate water-year runoff
windowSize = 12;
b = ones(1,windowSize);
a = 1;
Rs_S0_12 = filter(b, a, Rs_S0, [], 1);
Rs_S0_12(1:(windowSize-1), :, :) = NaN;
Rs_S1_12 = filter(b, a, Rs_S1, [], 1);
Rs_S1_12(1:(windowSize-1), :, :) = NaN;
Rs_S2_12 = filter(b, a, Rs_S2, [], 1);
Rs_S2_12(1:(windowSize-1), :, :) = NaN;
Rs_S3_12 = filter(b, a, Rs_S3, [], 1);
Rs_S3_12(1:(windowSize-1), :, :) = NaN;

Rs_S0_WY = Rs_S0_12(mos==9, :, :); 
Rs_S1_WY = Rs_S1_12(mos==9, :, :); 
Rs_S2_WY = Rs_S2_12(mos==9, :, :); 
Rs_S3_WY = Rs_S3_12(mos==9, :, :); 

% Aggregate to basin scale runoff
Q_S0_WY = sum((0.001 .* Rs_S0_WY(:, MRBidx) .* repmat(area(MRBidx)', size(Rs_S0_WY, 1), 1)) / (365 * 24 * 60 * 60), 2);
Q_S1_WY = sum((0.001 .* Rs_S1_WY(:, MRBidx) .* repmat(area(MRBidx)', size(Rs_S1_WY, 1), 1)) / (365 * 24 * 60 * 60), 2);
Q_S2_WY = sum((0.001 .* Rs_S2_WY(:, MRBidx) .* repmat(area(MRBidx)', size(Rs_S2_WY, 1), 1)) / (365 * 24 * 60 * 60), 2);
Q_S3_WY = sum((0.001 .* Rs_S3_WY(:, MRBidx) .* repmat(area(MRBidx)', size(Rs_S3_WY, 1), 1)) / (365 * 24 * 60 * 60), 2);

clear b a windowSize WB12 P12 E12;

%% Read in flow
gages = {'BNMO','FTRA','HEMO','MKC','NCNE','SUX','BIGSIOUX','GASCONADE','JAMES','LOSAGE','WHITE'};
T = readtable('../data/LMR_gages_WY_mean_flow.xlsx', ...
    'ReadVariableNames',false, 'Range','D2:N122');
T.Properties.VariableNames = gages;
T.Year = [1898:2018]';

%% Plot model vs. measurement comparison
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6.5 2.];
ylim = [0 8000];

subplot(1,5,1:3)
plot(T.Year, T.HEMO, 'k-', 'LineWidth',1.2)
hold on;
plot(year, Q_S3_WY, '-', 'LineWidth',1.2, 'Color',[0.5 0.5 0.5])
ax1 = gca;
set(ax1, 'XLim',[1900 2021],'YLim',ylim,'TickDir','out')
box off;
ylabel('Mean Q (m^{3} s^{-1})')
ax1.Position(1) = 0.1;
ax1.Position(2) = 0.22;
ax1.Position(3) = 0.46;
ax1.Position(4) = 0.76;
lgd = legend(ax1, 'Measured','Modeled', 'Location','northwest');
lgd.Position(2) = 0.84;

subplot(1,5,4:5)
[~,ia,ib] = intersect(year, T.Year);
plot(T.HEMO(ib), Q_S3_WY(ia), 'ko')
hold on;
plot(ylim, ylim, 'k-')
ax2 = gca;
set(ax2, 'XLim',ylim, 'YLim',ylim, 'TickDir','out')
box off;
ax2.Position(1) = 0.68;
ax2.Position(2) = 0.22;
ax2.Position(4) = 0.76;
xlabel('Measured mean Q (m^{-3} s^{-1})')
ylabel('Modeled mean Q (m^{-3} s^{-1})')
r2 = corr(T.HEMO(ib), Q_S3_WY(ia), 'rows','pairwise') ^ 2;
text(ylim(1)+diff(ylim)*0.05, ylim(2), sprintf('r^{2} = %.02f', r2), 'VerticalAlignment','top')

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','../output/wbm-hemo-Q.tif')
close all;

% %% bias correct
% [~,ia,ib] = intersect(year, T.Year);
% x = Q_S3_WY(ia); y = T.HEMO(ib);
% mx = mean(x(~isnan(x) & ~isnan(y)));
% my = mean(y(~isnan(x) & ~isnan(y)));
% sx = std(x(~isnan(x) & ~isnan(y)));
% sy = std(y(~isnan(x) & ~isnan(y)));
% Q_S0_WY_sc = my + sy*((Q_S0_WY - mx)/sx);
% Q_S1_WY_sc = my + sy*((Q_S1_WY - mx)/sx);
% Q_S2_WY_sc = my + sy*((Q_S2_WY - mx)/sx);
% Q_S3_WY_sc = my + sy*((Q_S3_WY - mx)/sx);

%% save data
cd ..;
save('./data/McCabeWilliams_WaterBudget.mat', 'year', 'lat', 'lon', 'MRBidx',...
    'Q_S0_WY', 'Q_S1_WY', 'Q_S2_WY', 'Q_S3_WY', 'Rs_S0_WY', "Rs_S1_WY", "Rs_S2_WY", "Rs_S3_WY", '-v7.3');

%% plot contributions
plot(1901:2021, Q_S3_WY - Q_S2_WY, 'b-')
hold on;
plot(1901:2021, Q_S2_WY - Q_S1_WY, 'r-')
plot(1901:2021, Q_S1_WY - Q_S0_WY, 'k-')
plot(year, Q_S3_WY - Q_S0_WY, '-', 'Color',[0.4 0.4 0.4], 'LineWidth',1.5)


