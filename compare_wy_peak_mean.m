Twy = readtable('./data/Missouri_Big6_UpdatedUnregulatedFlow_DailyMonthlyWY_avgCMS.xlsx',...
    'Sheet','WaterYear_cms');
Twy.Properties.VariableNames = {'Year','FTRA','SUX','NCNE','MKC','MNMO','HEMO'};
Twy.HEMO_peak = NaN(size(Twy,1),1);

Tdy = readtable('./data/Missouri_Big6_UpdatedUnregulatedFlow_DailyMonthlyWY_avgCMS.xlsx',...
    'Sheet','Daily_cms');
Tdy.Properties.VariableNames = {'Year','Month','Day','FTRA','SUX','NCNE','MKC','MNMO','HEMO'};

% Get water-year peak flows
yrs = Twy.Year;
for i = 1:length(yrs)
    idx = find(Tdy.Year==(yrs(i)-1) & Tdy.Month==10 & Tdy.Day==1):find(Tdy.Year==yrs(i) & Tdy.Month==9 & Tdy.Day==30);
    Twy.HEMO_peak(i) = max(Tdy.HEMO(idx));
end

% Make figure
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 4 3.5];

s1 = scatter(Twy.HEMO, Twy.HEMO_peak, 30, 'k','filled');
ax = ancestor(s1, 'axes');
ax.YAxis.Exponent = 0;
xlabel('Mean water-year flow (m^{3} s^{-1})')
ylabel('Peak water-year flow (m^{3} s^{-1})')
set(gca, 'TickDir','out')
[r,p] = corr(Twy.HEMO, Twy.HEMO_peak);
xlim = ax.XLim;
ylim = ax.YLim;
text(xlim(1)+0.02*diff(xlim), ylim(2),['{\itR} =  ',sprintf('%.2f', r)]);
text(xlim(1)+0.02*diff(xlim), ylim(2)-0.05*diff(ylim), '{\itp} < 0.001')

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r600','./output/wy-mean-vs-peak.tif')
close all;
