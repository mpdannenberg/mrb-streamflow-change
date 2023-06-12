%% Basic script parameters
syear = 1901;
eyear = 2010;
load ./data/MRB_subregions_GP;
latlim = [32 50];
lonlim = [-115 -85];
windowSize = 12;
b = ones(1,windowSize);
a = 1;

%% Read in flow
gages = {'BNMO','FTRA','HEMO','MKC','NCNE','SUX','BIGSIOUX','GASCONADE','JAMES','LOSAGE','WHITE'};
T = readtable('./data/LMR_gages_WY_mean_flow.xlsx', ...
    'ReadVariableNames',false, 'Range','D2:N122');
T.Properties.VariableNames = gages;
T.Year = [1898:2018]';

%% Calculate water-year climatic water balance (CWB) and aridity index (AI)
clatlim = [32 50];
clonlim = [-115 -85];
yr = reshape(repmat(1901:2019, 12, 1), [], 1);
mos = repmat(1:12, 1, length(1901:2019))';
lat = ncread("data\cru_ts4.04.1901.2019.pre.dat.nc",'lat'); latidx = lat >= clatlim(1) & lat <= clatlim(2);
lon = ncread("data\cru_ts4.04.1901.2019.pre.dat.nc",'lon'); lonidx = lon >= clonlim(1) & lon <= clonlim(2);
P = permute(ncread("data\cru_ts4.04.1901.2019.pre.dat.nc",'pre'), [3 2 1]); P = P(:, latidx, lonidx);

P12 = filter(b, a, P, [], 1);
P12(1:(windowSize-1), :, :) = NaN;
P = P12(mos==9 & yr>=syear & yr<=eyear, :, :) * (1/1000); % convert mm year-1 to m year-1

[nt, ny, nx] = size(P);
clear P12 yr mos lat lon;

%% Get MsTMIP model info
cd('D:/Data_Analysis/MsTMIP/')
models = {'CLM4','CLM4VIC','DLEM','GTEC','ISAM','LPJ-wsl',...
    'ORCHIDEE-LSCE','SiBCASA','VEGAS2.1','VISIT'}; % exclude CLASS-CTEM-N due to negative water balance (almost 100% of the time) and SiB3 due to negative water balance (most of the time) in SG2 and SG3 scenarios
nk = length(models);

info = ncinfo([models{1},'_BG1_Monthly_Evap.nc4']);
lat = double(ncread([models{1},'_BG1_Monthly_Evap.nc4'], 'lat')); latidx = lat >= latlim(1) & lat <= latlim(2); 
lon = double(ncread([models{1},'_BG1_Monthly_Evap.nc4'], 'lon')); lonidx = lon >= lonlim(1) & lon <= lonlim(2); 
yr = reshape(repmat(syear:eyear, 12, 1), [], 1);
mo = reshape(repmat(1:12, nt, 1)', [], 1);
eom = eomday(yr, mo); % number of days per month

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

clear SR* LAT LON LatLon IN ON idx e;

%% Loop through models and get total ET
ETtot = NaN(nt, nk);
ET_RG1 = NaN(nt, ny, nx, nk);
ET_SG1 = NaN(nt, ny, nx, nk);
ET_SG2 = NaN(nt, ny, nx, nk);
ET_SG3 = NaN(nt, ny, nx, nk);

for k = 1:nk
    
    % RG1: baseline simulation
    evap = permute(ncread([models{k},'_RG1_Monthly_Evap.nc4'], 'Evap'), [3 2 1]);
    evap = evap(:, latidx, lonidx) .* (repmat(eom, 1, ny, nx)*24*60*60) * 0.001; % subset and convert from kg m-2 s-1 --> m3 m-2 month-1
    ET = filter(b, a, evap, [], 1); % 12-month sums (m3 m-2 month-1 --> m3 m-2 yr-1 [or m year-1])
    ET(1:(windowSize-1), :, :) = NaN;
    ET_RG1(:,:,:,k) = ET(mo==9,:,:);
    
    % SG1: climate-only simulation
    evap = permute(ncread([models{k},'_SG1_Monthly_Evap.nc4'], 'Evap'), [3 2 1]);
    evap = evap(:, latidx, lonidx) .* (repmat(eom, 1, ny, nx)*24*60*60) * 0.001; % subset and convert from kg m-2 s-1 --> m3 m-2 month-1
    ET = filter(b, a, evap, [], 1); % 12-month sums (m3 m-2 month-1 --> m3 m-2 yr-1 [or m year-1])
    ET(1:(windowSize-1), :, :) = NaN;
    ET_SG1(:,:,:,k) = ET(mo==9,:,:);
    
    % SG2: climate + LULCC
    evap = permute(ncread([models{k},'_SG2_Monthly_Evap.nc4'], 'Evap'), [3 2 1]);
    evap = evap(:, latidx, lonidx) .* (repmat(eom, 1, ny, nx)*24*60*60) * 0.001; % subset and convert from kg m-2 s-1 --> m3 m-2 month-1
    ET = filter(b, a, evap, [], 1); % 12-month sums (m3 m-2 month-1 --> m3 m-2 yr-1 [or m year-1])
    ET(1:(windowSize-1), :, :) = NaN;
    ET_SG2(:,:,:,k) = ET(mo==9,:,:);
    
    % SG3: climate + LULCC + CO2
    evap = permute(ncread([models{k},'_SG3_Monthly_Evap.nc4'], 'Evap'), [3 2 1]);
    evap = evap(:, latidx, lonidx) .* (repmat(eom, 1, ny, nx)*24*60*60) * 0.001; % subset and convert from kg m-2 s-1 --> m3 m-2 month-1
    ET = filter(b, a, evap, [], 1); % 12-month sums (m3 m-2 month-1 --> m3 m-2 yr-1 [or m year-1])
    ET(1:(windowSize-1), :, :) = NaN;
    ET_SG3(:,:,:,k) = ET(mo==9,:,:);
    
end

%% Calculate water balance using CRU precip and MsTMIP ET
cd('D:\Publications\NSF_LMRB_CFR')
Ptot = sum(repmat(area(MRBidx)', nt, 1) .* P(:,MRBidx), 2);
Pmean = sum(repmat(area(MRBidx)', nt, 1) .* repmat(mean(P(1:30,MRBidx),'omitnan'),nt,1), 2);
ETtot_RG1 = NaN(nt, nk);
for k = 1:nk
    temp = ET_RG1(:,:,:,k);
    ETtot_RG1(:,k) = sum(repmat(area(MRBidx)', nt, 1) .* temp(:,MRBidx), 2);
end
ETtot_SG1 = NaN(nt, nk);
for k = 1:nk
    temp = ET_SG1(:,:,:,k);
    ETtot_SG1(:,k) = sum(repmat(area(MRBidx)', nt, 1) .* temp(:,MRBidx), 2);
end
ETtot_SG2 = NaN(nt, nk);
for k = 1:nk
    temp = ET_SG2(:,:,:,k);
    ETtot_SG2(:,k) = sum(repmat(area(MRBidx)', nt, 1) .* temp(:,MRBidx), 2);
end
ETtot_SG3 = NaN(nt, nk);
for k = 1:nk
    temp = ET_SG3(:,:,:,k);
    ETtot_SG3(:,k) = sum(repmat(area(MRBidx)', nt, 1) .* temp(:,MRBidx), 2);
end

% convert from total m3 yr-1 to m3 s-1
WB_RG1 = (repmat(Pmean,1,nk) - ETtot_RG1) * (1/365) * (1/24) * (1/60) * (1/60);
WB_SG1 = (repmat(Ptot,1,nk) - ETtot_SG1) * (1/365) * (1/24) * (1/60) * (1/60);
WB_SG2 = (repmat(Ptot,1,nk) - ETtot_SG2) * (1/365) * (1/24) * (1/60) * (1/60);
WB_SG3 = (repmat(Ptot,1,nk) - ETtot_SG3) * (1/365) * (1/24) * (1/60) * (1/60);

% %% scale to match mean/variance of observed streamflow
% WB_RG1_sc = NaN(size(WB_RG1));
% WB_SG1_sc = NaN(size(WB_SG1));
% WB_SG2_sc = NaN(size(WB_SG2));
% WB_SG3_sc = NaN(size(WB_SG3));
% for k = 1:nk
%     wb = WB_SG3(:,k);
%     [~,ia,ib] = intersect(syear:eyear, T.Year);
%     
%     x = wb(ia); y = T.HEMO(ib);
%     mx = mean(x(~isnan(x) & ~isnan(y)));
%     my = mean(y(~isnan(x) & ~isnan(y)));
%     sx = std(x(~isnan(x) & ~isnan(y)));
%     sy = std(y(~isnan(x) & ~isnan(y)));
% 
%     WB_RG1_sc(:,k) = my + sy*((WB_RG1(:,k) - mx)/sx);
%     WB_SG1_sc(:,k) = my + sy*((WB_SG1(:,k) - mx)/sx);
%     WB_SG2_sc(:,k) = my + sy*((WB_SG2(:,k) - mx)/sx);
%     WB_SG3_sc(:,k) = my + sy*((WB_SG3(:,k) - mx)/sx);
%     
% %     mdl = fitlm(wb(ia), T.HEMO(ib));
% %     mdl.predict
% end
% WB_RG1 = WB_RG1_sc;
% WB_SG1 = WB_SG1_sc;
% WB_SG2 = WB_SG2_sc;
% WB_SG3 = WB_SG3_sc;

%% Save data
year = syear:eyear;
lat = lat(latidx);
lon = lon(lonidx);
save('./data/MsTMIP_WaterBudget.mat', 'year', 'lat', 'lon', 'models',...
    'ET_RG1','ET_SG1','ET_SG2','ET_SG3','WB_RG1','WB_SG1','WB_SG2','WB_SG3', 'P', '-v7.3');

% %% make figures
% h = figure('Color','w');
% h.Units = 'inches';
% h.Position = [1 1 5 7];
% 
% subplot(3,1,1)
% plot(repmat([syear:eyear]', 1, nk), WB_SG1)
% set(gca, 'XTickLabel','', 'TickDir','out', 'XLim',[syear eyear], 'YLim',[0 15000])
% box off;
% hold on;
% plot(T.Year, T.HEMO, 'k-', 'LineWidth',1.5)
% title('SG1')
% 
% subplot(3,1,2)
% plot(repmat([syear:eyear]', 1, nk), WB_SG2)
% set(gca, 'XTickLabel','', 'TickDir','out', 'XLim',[syear eyear], 'YLim',[0 15000])
% box off;
% hold on;
% plot(T.Year, T.HEMO, 'k-', 'LineWidth',1.5)
% title('SG2')
% ylabel('Streamflow (m^{3} s^{-1})', 'FontSize',11)
% 
% subplot(3,1,3)
% plot(repmat([syear:eyear]', 1, nk), WB_SG3)
% set(gca, 'TickDir','out', 'XLim',[syear eyear], 'YLim',[0 15000])
% box off;
% hold on;
% plot(T.Year, T.HEMO, 'k-', 'LineWidth',1.5)
% title('SG3')
% 
% %% make figure showing unscaled and rescaled WB vs. HEMO flow
% h = figure('Color','w');
% h.Units = 'inches';
% h.Position = [1 1 5 6];
% 
% subplot(2,1,1)
% plot(syear:eyear, WB_SG3)
% set(gca, 'XTickLabel','', 'TickDir','out', 'XLim',[syear eyear], 'YLim',[0 15000])
% box off;
% hold on;
% plot(T.Year, T.HEMO, 'k-', 'LineWidth',1.5)
% title('unscaled')
% hold off;
% ylabel('Streamflow (m^{3} s^{-1})', 'FontSize',11)
% 
% subplot(2,1,2)
% plot(syear:eyear, WB_SG3_sc)
% set(gca, 'XTickLabel','', 'TickDir','out', 'XLim',[syear eyear], 'YLim',[0 6000])
% box off;
% hold on;
% plot(T.Year, T.HEMO, 'k-', 'LineWidth',1.5)
% title('scaled')
% ylabel('Streamflow (m^{3} s^{-1})', 'FontSize',11)
% 
% 
% %% plot change due to each factor
% h = figure('Color','w');
% h.Units = 'inches';
% h.Position = [1 1 5 7];
% 
% subplot(3,1,1)
% plot(repmat([syear:eyear]', 1, nk), WB_SG1 - WB_RG1)
% set(gca, 'XTickLabel','', 'TickDir','out', 'XLim',[syear eyear], 'YLim',[-4000 4000])
% box off;
% hold on;
% % plot(T.Year, T.HEMO, 'k-', 'LineWidth',1.5)
% title('Climate')
% plot(syear:eyear, median(WB_SG1 - WB_RG1, 2),'k-','LineWidth',1.5)
% 
% subplot(3,1,2)
% plot(repmat([syear:eyear]', 1, nk), WB_SG2 - WB_SG1)
% set(gca, 'XTickLabel','', 'TickDir','out', 'XLim',[syear eyear], 'YLim',[-4000 4000])
% box off;
% hold on;
% % plot(T.Year, T.HEMO, 'k-', 'LineWidth',1.5)
% title('LULCC')
% ylabel('Streamflow (m^{3} s^{-1})', 'FontSize',11)
% plot(syear:eyear, median(WB_SG2 - WB_SG1, 2),'k-','LineWidth',1.5)
% 
% subplot(3,1,3)
% plot(repmat([syear:eyear]', 1, nk), WB_SG3 - WB_SG2)
% set(gca, 'TickDir','out', 'XLim',[syear eyear], 'YLim',[-4000 4000])
% box off;
% hold on;
% % plot(T.Year, T.HEMO, 'k-', 'LineWidth',1.5)
% title('CO2')
% plot(syear:eyear, median(WB_SG3 - WB_SG2, 2),'k-','LineWidth',1.5)

