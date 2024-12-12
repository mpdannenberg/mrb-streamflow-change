
% basic parameters
syear = 1950;
eyear = 2021;

% filter parameters 
windowSize = 3;
a = 1;
b = ones(1,windowSize) / windowSize;

% full table setup
T = table('Size',[length(syear:eyear)*12 7], ...
    'VariableTypes',repmat({'double'},1,7), ...
    'VariableNames',{'Year','Month','AMO','NAO','Nino34','PDO','PNA'});
T.Year = reshape(repmat(syear:eyear, 12, 1), [], 1);
T.Month = repmat([1:12]', length(syear:eyear), 1);

% load indices
amo = readtable('./data/teleconnecctions/amo_monthly.txt', 'TreatAsMissing',{'-999'}); % https://climatedataguide.ucar.edu/sites/default/files/2022-03/amo_monthly.txt
pna = readtable('./data/teleconnecctions/norm.pna.monthly.txt'); % https://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/norm.pna.monthly.b5001.current.ascii.table
nino34 = readtable('./data/teleconnecctions/oni.txt'); % https://psl.noaa.gov/data/correlation/oni.data
nao = readtable('./data/teleconnecctions/nao_pc_monthly.txt'); % https://climatedataguide.ucar.edu/sites/default/files/2024-04/nao_pc_monthly.txt

% add to table
T.AMO = filter(b, a, reshape(amo{:,2:end}', [], 1));
T.NAO = filter(b, a, reshape(nao{:,2:end}', [], 1));
T.Nino34 = reshape(nao{:,2:end}', [], 1); % Nino3.4 already 3-month average (ONI)
T.PNA = filter(b, a, reshape(pna{:,2:end}', [], 1));

% MW11 water budget
load('./data/McCabeWilliams_WaterBudget.mat')
Q_nat = Q_S0_WY(year >= syear & year <= eyear);
clear year lat lon MRBidx Q_S* Rs_*;

% Correlation "natural" water budget with teleconnection indices
[r(1,:), p(1,:)] = corr(Q_nat, lagmatrix(reshape(T.AMO,12,[])', [1 0]), "rows","pairwise");
[r(end+1,:), p(end+1,:)] = corr(Q_nat, lagmatrix(reshape(T.NAO,12,[])', [1 0]), "rows","pairwise");
[r(end+1,:), p(end+1,:)] = corr(Q_nat, lagmatrix(reshape(T.Nino34,12,[])', [1 0]), "rows","pairwise");
[r(end+1,:), p(end+1,:)] = corr(Q_nat, lagmatrix(reshape(T.PNA,12,[])', [1 0]), "rows","pairwise");
[~,pfdr] = fdr(reshape(p(:,1:21),1,[]), 0.1);

% make figure
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 4.5 5.5];
ax = tight_subplot(4,1,0.04,0.05,[0.12 0.04]);
mlabels = {'j','f','m','a','m','j','j','a','s','o','n','d','J','F','M','A','M','J','J','A','S'};

axes(ax(1))
b = bar(r(1,1:21), 'FaceColor','flat');
b.EdgeColor = 'none';
b.CData(p(1,1:21) > pfdr, :) = repmat([0.8 0.8 0.8],sum(p(1,1:21) > pfdr),1);
b.CData(p(1,1:21) <= pfdr, :) = repmat([0.3 0.3 0.3],sum(p(1,1:21) <= pfdr),1);
set(gca, 'XLim',[0.5 21.5], 'TickDir','out', 'XColor','none', 'YColor','k', 'YGrid','on')
box off;
ylabel('{\itR}(Q_{Nat}, AMO)')
ylim = get(gca, 'YLim');
text(1,ylim(1)+0.05*diff(ylim),'a','FontSize',12,'FontWeight','bold','VerticalAlignment','baseline')

axes(ax(2))
b = bar(r(2,1:21), 'FaceColor','flat');
b.EdgeColor = 'none';
b.CData(p(2,1:21) > pfdr, :) = repmat([0.8 0.8 0.8],sum(p(2,1:21) > pfdr),1);
b.CData(p(2,1:21) <= pfdr, :) = repmat([0.3 0.3 0.3],sum(p(2,1:21) <= pfdr),1);
set(gca, 'XLim',[0.5 21.5], 'TickDir','out', 'XColor','none', 'YColor','k', 'YGrid','on')
box off;
ylabel('{\itR}(Q_{Nat}, NAO)')
ylim = get(gca, 'YLim');
text(1,ylim(1)+0.05*diff(ylim),'b','FontSize',12,'FontWeight','bold','VerticalAlignment','baseline')

axes(ax(3))
b = bar(r(3,1:21), 'FaceColor','flat');
b.EdgeColor = 'none';
b.CData(p(3,1:21) > pfdr, :) = repmat([0.8 0.8 0.8],sum(p(3,1:21) > pfdr),1);
b.CData(p(3,1:21) <= pfdr, :) = repmat([0.3 0.3 0.3],sum(p(3,1:21) <= pfdr),1);
set(gca, 'XLim',[0.5 21.5], 'TickDir','out', 'XColor','none', 'YColor','k', 'YGrid','on')
box off;
ylabel('{\itR}(Q_{Nat}, Niño3.4)')
ylim = get(gca, 'YLim');
text(1,ylim(1)+0.05*diff(ylim),'c','FontSize',12,'FontWeight','bold','VerticalAlignment','baseline')

axes(ax(4))
b = bar(r(4,1:21), 'FaceColor','flat');
b.EdgeColor = 'none';
b.CData(p(4,1:21) > pfdr, :) = repmat([0.8 0.8 0.8],sum(p(4,1:21) > pfdr),1);
b.CData(p(4,1:21) <= pfdr, :) = repmat([0.3 0.3 0.3],sum(p(4,1:21) <= pfdr),1);
set(gca, 'XLim',[0.5 21.5], 'TickDir','out', 'XColor','none', 'YColor','k', 'YGrid','on')
box off;
ylim = get(gca, 'YLim');
for i=1:21; text(i,ylim(1),mlabels{i},'HorizontalAlignment','center', 'VerticalAlignment','top', 'FontSize',10); end
ylabel('{\itR}(Q_{Nat}, PNA)')
text(1,ylim(1)+0.05*diff(ylim),'d','FontSize',12,'FontWeight','bold','VerticalAlignment','baseline')

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r600','./output/qnat_teleconnection_correlations.tif')
close all;
