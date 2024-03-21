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

% contributions by subregion
Clim = permute(Q_SG1-Q_RG1, [1 3 2]);
LULCC = permute(Q_SG2-Q_SG1, [1 3 2]);
CO2 = permute(Q_SG3-Q_SG2, [1 3 2]);
All = permute(Q_SG3-Q_RG1, [1 3 2]);

T = table('Size',[7 5], 'VariableTypes',repmat({'string'},1,5),...
    'VariableNames',{'Region','dQ','dQ_clim','dQ_lulcc','dQ_co2'});
T.Region = {'MRB','1','2','3','4','5','6'}';

% All
m = squeeze(mean(All(idx, :, :), [1 2]));
for i=1:7
    r = corr(All(:,:,i),'rows','pairwise'); r = mean(r(~eye(size(r))).^2);
    s = squeeze(std(All(idx, :, i), 0, [1 2])) / sqrt(sum(idx)*nk*(1-r));
    T.dQ(i) = [sprintf('%0.1f ',m(i)),char(177),sprintf(' %0.1f',1.96*s)]; 
end

% Clim
m = squeeze(mean(Clim(idx, :, :), [1 2]));
for i=1:7
    r = corr(Clim(:,:,i),'rows','pairwise'); r = mean(r(~eye(size(r))).^2);
    s = squeeze(std(Clim(idx, :, i), 0, [1 2])) / sqrt(sum(idx)*nk*(1-r));
    T.dQ_clim(i) = [sprintf('%0.1f ',m(i)),char(177),sprintf(' %0.1f',1.96*s)]; 
end

% LULCC
m = squeeze(mean(LULCC(idx, :, :), [1 2]));
for i=1:7
    r = corr(LULCC(:,:,i),'rows','pairwise'); r = mean(r(~eye(size(r))).^2);
    s = squeeze(std(LULCC(idx, :, i), 0, [1 2])) / sqrt(sum(idx)*nk*(1-r));
    T.dQ_lulcc(i) = [sprintf('%0.1f ',m(i)),char(177),sprintf(' %0.1f',1.96*s)]; 
end

% CO2
m = squeeze(mean(CO2(idx, :, :), [1 2]));
for i=1:7
    r = corr(CO2(:,:,i),'rows','pairwise'); r = mean(r(~eye(size(r))).^2);
    s = squeeze(std(CO2(idx, :, i), 0, [1 2])) / sqrt(sum(idx)*nk*(1-r));
    T.dQ_co2(i) = [sprintf('%0.1f ',m(i)),char(177),sprintf(' %0.1f',1.96*s)]; 
end

writetable(T, './output/table-s5-mstmip-effects.xlsx');

