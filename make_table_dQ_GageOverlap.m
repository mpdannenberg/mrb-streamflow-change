% TWO ISSUES HERE (and fix in test_variance... file as well)
    % 1) square R before averaging? revisit...
    % 2) multiple by number of models to get Ns in MsTMIP section


dQ = table('Size',[7 2],'VariableTypes',{'string','string'},...
    'VariableNames',{'Var','dQ'});
dQens = table('Size',[10 4],'VariableTypes',{'string','string','string','string'},...
    'VariableNames',{'Model','dQ_clim','dQ_lulcc','dQ_co2'});

clr = wesanderson('fantasticfox1');
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6.5 3];

% gage
T = readtable('./data/Missouri_Big6_UpdatedUnregulatedFlow_DailyMonthlyWY_avgCMS.xlsx',...
    'Sheet','WaterYear_cms');
T.Properties.VariableNames = {'Year','FTRA','SUX','NCNE','MKC','MNMO','HEMO'};
y0 = T.HEMO(T.Year>=1931 & T.Year<=1960);
y1 = T.HEMO(T.Year>=1981 & T.Year<=2010);
x0 = mean(y0); s0 = std(y0); n0 = length(y0);
x1 = mean(y1); s1 = std(y1); n1 = length(y1);
sp = sqrt( ((n0-1)*s0^2 + (n1-1)*s1^2) / (n0+n1-2) );
se = sp * sqrt(1/n0 + 1/n1);
plot(T.Year, T.HEMO-x0, 'k-', 'LineWidth',1.5)
hold on;
plot([1930 2010],[0 0],'k-','LineWidth',0.5)
% plot([1931 2010], [x0 x0], 'k-', 'LineWidth',0.5)
plot([1981 2010], [x1 x1]-x0, 'k-', 'LineWidth',0.75)
set(gca, 'TickDir','out', 'XLim', [1930 2010])
box off;
ylabel('Q (m^{3} s^{-1})')
ax = gca;
ax.Position(1) = 0.11;
ax.Position(3) = 0.73;
text(2010.5, x1-x0, ['+',sprintf('%0.1f',x1-x0),' m^{3} s^{-1}'], 'VerticalAlignment','middle')
dQ.Var(1) = 'dQ_gage';
dQ.dQ(1) = [sprintf('%0.1f ',x1-x0),char(177),sprintf(' %0.1f',1.96*se)];

% WBM
load('./data/McCabeWilliams_WaterBudget.mat')
y0 = Q_S0_WY(year>=1931 & year<=1960);
y1 = Q_S0_WY(year>=1981 & year<=2010);
x0 = mean(y0); s0 = std(y0); n0 = length(y0);
x1 = mean(y1); s1 = std(y1); n1 = length(y1);
sp = sqrt( ((n0-1)*s0^2 + (n1-1)*s1^2) / (n0+n1-2) );
se = sp * sqrt(1/n0 + 1/n1);
plot(year, Q_S0_WY-x0, '-', 'Color',clr(2,:), 'LineWidth',1.5)
plot([1981 2010], [x1 x1]-x0, '-', 'Color',clr(2,:), 'LineWidth',0.75)
% text(2010, mean([x0 x1]), ['+',sprintf('%0.1f',x1-x0),' m^{3} s^{-1}'], 'Color',clr(2,:))
dQ.Var(2) = 'dQ_climate_nat';
dQ.dQ(2) = [sprintf('%0.1f ',x1-x0),char(177),sprintf(' %0.1f',1.96*se)];

y0 = Q_S3_WY(year>=1931 & year<=1960) - Q_S0_WY(year>=1931 & year<=1960);
y1 = Q_S3_WY(year>=1981 & year<=2010) - Q_S0_WY(year>=1981 & year<=2010);
x0 = mean(y0); s0 = std(y0); n0 = length(y0);
x1 = mean(y1); s1 = std(y1); n1 = length(y1);
sp = sqrt( ((n0-1)*s0^2 + (n1-1)*s1^2) / (n0+n1-2) );
se = sp * sqrt(1/n0 + 1/n1);
plot(year, Q_S3_WY-Q_S0_WY-x0, '--', 'Color',clr(2,:), 'LineWidth',1.5)
dQ.Var(3) = 'dQ_climate_anth';
dQ.dQ(3) = [sprintf('%0.1f ',x1-x0),char(177),sprintf(' %0.1f',1.96*se)];

y0 = Q_S3_WY(year>=1931 & year<=1960);
y1 = Q_S3_WY(year>=1981 & year<=2010);
x0 = mean(y0); s0 = std(y0); n0 = length(y0);
x1 = mean(y1); s1 = std(y1); n1 = length(y1);
sp = sqrt( ((n0-1)*s0^2 + (n1-1)*s1^2) / (n0+n1-2) );
se = sp * sqrt(1/n0 + 1/n1);
dQ.Var(4) = 'dQ_climate_wbm';
dQ.dQ(4) = [sprintf('%0.1f ',x1-x0),char(177),sprintf(' %0.1f',1.96*se)];

% MsTMIP
% Calculating an effective sample size for the number of models by:
% Neff = N / (1 + r), where r is the mean between-model correlation
% via ChatGPT, who cited Heeringa et al. (2017), Applied Survey Data
% Analysis, CRC Press. [probably worth checking that]

% Actually, no. Using variance deflation factor (1-r^2) from Clifton et al.
% 2019. Could also do Dannenberg et al. 2015? But way too conservative?

% Actually, double no. Do ensemble mean and each one individually
load('./output/MsTMIP-contributions.mat')
dQens.Model = models'; nk = length(models);
Clim = Q_SG1 - Q_RG1; % r_clim = corr(Clim,'rows','pairwise'); r_clim = mean(r_clim(~eye(size(r_clim))));
LULCC = Q_SG2 - Q_SG1; % r_lulcc = corr(LULCC,'rows','pairwise');  r_lulcc = mean(r_lulcc(~eye(size(r_lulcc))));
CO2 = Q_SG3 - Q_SG2; % r_co2 = corr(CO2,'rows','pairwise');  r_co2 = mean(r_co2(~eye(size(r_co2))));

% y0 = reshape(Clim(year>=1931 & year<=1960,:),[],1);
% y1 = reshape(Clim(year>=1981 & year<=2010,:),[],1);
y0 = mean(Clim(year>=1931 & year<=1960,:), 2);
y1 = mean(Clim(year>=1981 & year<=2010,:), 2);
x0 = mean(y0); s0 = std(y0); n0 = length(y0); % * (1 - r_clim^2);
x1 = mean(y1); s1 = std(y1); n1 = length(y1); % * (1 - r_clim^2);
sp = sqrt( ((n0-1)*s0^2 + (n1-1)*s1^2) / (n0+n1-2) );
se = sp * sqrt(1/n0 + 1/n1);
plot(year, mean(Clim,2)-x0, '-', 'Color',clr(1,:), 'LineWidth',1.5)
plot([1981 2010], [x1 x1]-x0, '-', 'Color',clr(1,:), 'LineWidth',0.75)
dQ.Var(5) = 'dQ_climate_mstmip';
dQ.dQ(5) = [sprintf('%0.1f ',x1-x0),char(177),sprintf(' %0.1f',1.96*se)];
for i = 1:nk
    y0 = Clim(year>=1931 & year<=1960,i);
    y1 = Clim(year>=1981 & year<=2010,i);
    x0 = mean(y0); s0 = std(y0); n0 = length(y0); % * (1 - r_clim^2);
    x1 = mean(y1); s1 = std(y1); n1 = length(y1); % * (1 - r_clim^2);
    sp = sqrt( ((n0-1)*s0^2 + (n1-1)*s1^2) / (n0+n1-2) );
    se = sp * sqrt(1/n0 + 1/n1);
    dQens.dQ_clim(i) = [sprintf('%0.1f ',x1-x0),char(177),sprintf(' %0.1f',1.96*se)];
end

% y0 = reshape(LULCC(year>=1931 & year<=1960,:),[],1);
% y1 = reshape(LULCC(year>=1981 & year<=2010,:),[],1);
y0 = mean(LULCC(year>=1931 & year<=1960,:), 2);
y1 = mean(LULCC(year>=1981 & year<=2010,:), 2);
x0 = mean(y0); s0 = std(y0); n0 = length(y0); % * (1 - r_lulcc^2);
x1 = mean(y1); s1 = std(y1); n1 = length(y1); % * (1 - r_lulcc^2);
sp = sqrt( ((n0-1)*s0^2 + (n1-1)*s1^2) / (n0+n1-2) );
se = sp * sqrt(1/n0 + 1/n1);
plot(year, mean(LULCC,2)-x0, '-', 'Color',clr(3,:), 'LineWidth',1.5)
plot([1981 2010], [x1 x1]-x0, '-', 'Color',clr(3,:), 'LineWidth',0.75)
dQ.Var(6) = 'dQ_lulcc_mstmip';
dQ.dQ(6) = [sprintf('%0.1f ',x1-x0),char(177),sprintf(' %0.1f',1.96*se)];
for i = 1:nk
    y0 = LULCC(year>=1931 & year<=1960,i);
    y1 = LULCC(year>=1981 & year<=2010,i);
    x0 = mean(y0); s0 = std(y0); n0 = length(y0); % * (1 - r_clim^2);
    x1 = mean(y1); s1 = std(y1); n1 = length(y1); % * (1 - r_clim^2);
    sp = sqrt( ((n0-1)*s0^2 + (n1-1)*s1^2) / (n0+n1-2) );
    se = sp * sqrt(1/n0 + 1/n1);
    dQens.dQ_lulcc(i) = [sprintf('%0.1f ',x1-x0),char(177),sprintf(' %0.1f',1.96*se)];
end

% y0 = reshape(CO2(year>=1931 & year<=1960,:),[],1);
% y1 = reshape(CO2(year>=1981 & year<=2010,:),[],1);
y0 = mean(CO2(year>=1931 & year<=1960,:), 2);
y1 = mean(CO2(year>=1981 & year<=2010,:), 2);
x0 = mean(y0); s0 = std(y0); n0 = length(y0); % * (1 - r_co2^2);
x1 = mean(y1); s1 = std(y1); n1 = length(y1); % * (1 - r_co2^2);
sp = sqrt( ((n0-1)*s0^2 + (n1-1)*s1^2) / (n0+n1-2) );
se = sp * sqrt(1/n0 + 1/n1);
plot(year, mean(CO2,2)-x0, '-', 'Color',clr(4,:), 'LineWidth',1.5)
plot([1981 2010], [x1 x1]-x0, '-', 'Color',clr(4,:), 'LineWidth',0.75)
dQ.Var(7) = 'dQ_co2_mstmip';
dQ.dQ(7) = [sprintf('%0.1f ',x1-x0),char(177),sprintf(' %0.1f',1.96*se)];
for i = 1:nk
    y0 = CO2(year>=1931 & year<=1960,i);
    y1 = CO2(year>=1981 & year<=2010,i);
    x0 = mean(y0); s0 = std(y0); n0 = length(y0); % * (1 - r_clim^2);
    x1 = mean(y1); s1 = std(y1); n1 = length(y1); % * (1 - r_clim^2);
    sp = sqrt( ((n0-1)*s0^2 + (n1-1)*s1^2) / (n0+n1-2) );
    se = sp * sqrt(1/n0 + 1/n1);
    dQens.dQ_co2(i) = [sprintf('%0.1f ',x1-x0),char(177),sprintf(' %0.1f',1.96*se)];
end

writetable(dQ, './output/difference-in-means-instrumental.xlsx');
writetable(dQens, './output/difference-in-means-mstmip-models-instrumental.xlsx');
