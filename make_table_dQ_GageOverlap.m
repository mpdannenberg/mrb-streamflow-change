% Calculate change in mean and 95% CIs for each factor over instrumental
% period with gage overlap (1931-1960 vs. 1981-2010)
dQ = table('Size',[7 2],'VariableTypes',{'string','string'},...
    'VariableNames',{'Var','dQ'});
dQens = table('Size',[10 4],'VariableTypes',{'string','string','string','string'},...
    'VariableNames',{'Model','dQ_clim','dQ_lulcc','dQ_co2'});

Qbeta = table('Size',[7 2],'VariableTypes',{'string','string'},...
    'VariableNames',{'Var','Beta'});
Qbetaens = table('Size',[10 4],'VariableTypes',{'string','string','string','string'},...
    'VariableNames',{'Model','Beta_clim','Beta_lulcc','Beta_co2'});

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
dQ.Var(1) = 'dQ_gage';
dQ.dQ(1) = [sprintf('%0.1f ',x1-x0),char(177),sprintf(' %0.1f',1.96*se)];
mdl = fitlm(1931:2010, T.HEMO(T.Year>=1931 & T.Year<=2010));
Qbeta.Var(1) = 'Q_gage';
Qbeta.Beta(1) = [sprintf('%0.1f ',mdl.Coefficients.Estimate(2)),char(177),sprintf(' %0.1f',1.96*mdl.Coefficients.SE(2))];

% WBM
load('./data/McCabeWilliams_WaterBudget.mat')
y0 = Q_S0_WY(year>=1931 & year<=1960);
y1 = Q_S0_WY(year>=1981 & year<=2010);
x0 = mean(y0); s0 = std(y0); n0 = length(y0);
x1 = mean(y1); s1 = std(y1); n1 = length(y1);
sp = sqrt( ((n0-1)*s0^2 + (n1-1)*s1^2) / (n0+n1-2) );
se = sp * sqrt(1/n0 + 1/n1);
dQ.Var(2) = 'dQ_climate_nat';
dQ.dQ(2) = [sprintf('%0.1f ',x1-x0),char(177),sprintf(' %0.1f',1.96*se)];
y = Q_S0_WY;
mdl = fitlm(1931:2010, y(year>=1931 & year<=2010));
Qbeta.Var(2) = 'Q_climate_nat';
Qbeta.Beta(2) = [sprintf('%0.1f ',mdl.Coefficients.Estimate(2)),char(177),sprintf(' %0.1f',1.96*mdl.Coefficients.SE(2))];

y0 = Q_S3_WY(year>=1931 & year<=1960) - Q_S0_WY(year>=1931 & year<=1960);
y1 = Q_S3_WY(year>=1981 & year<=2010) - Q_S0_WY(year>=1981 & year<=2010);
x0 = mean(y0); s0 = std(y0); n0 = length(y0);
x1 = mean(y1); s1 = std(y1); n1 = length(y1);
sp = sqrt( ((n0-1)*s0^2 + (n1-1)*s1^2) / (n0+n1-2) );
se = sp * sqrt(1/n0 + 1/n1);
dQ.Var(3) = 'dQ_climate_anth';
dQ.dQ(3) = [sprintf('%0.1f ',x1-x0),char(177),sprintf(' %0.1f',1.96*se)];
y = Q_S3_WY - Q_S0_WY;
mdl = fitlm(1931:2010, y(year>=1931 & year<=2010));
Qbeta.Var(3) = 'Q_climate_anth';
Qbeta.Beta(3) = [sprintf('%0.1f ',mdl.Coefficients.Estimate(2)),char(177),sprintf(' %0.1f',1.96*mdl.Coefficients.SE(2))];

y0 = Q_S3_WY(year>=1931 & year<=1960);
y1 = Q_S3_WY(year>=1981 & year<=2010);
x0 = mean(y0); s0 = std(y0); n0 = length(y0);
x1 = mean(y1); s1 = std(y1); n1 = length(y1);
sp = sqrt( ((n0-1)*s0^2 + (n1-1)*s1^2) / (n0+n1-2) );
se = sp * sqrt(1/n0 + 1/n1);
dQ.Var(4) = 'dQ_climate_wbm';
dQ.dQ(4) = [sprintf('%0.1f ',x1-x0),char(177),sprintf(' %0.1f',1.96*se)];
y = Q_S3_WY;
mdl = fitlm(1931:2010, y(year>=1931 & year<=2010));
Qbeta.Var(4) = 'Q_climate_wbm';
Qbeta.Beta(4) = [sprintf('%0.1f ',mdl.Coefficients.Estimate(2)),char(177),sprintf(' %0.1f',1.96*mdl.Coefficients.SE(2))];

% MsTMIP
load('./output/MsTMIP-contributions.mat')
dQens.Model = models'; Qbetaens.Model = models'; nk = length(models);
Clim = Q_SG1 - Q_RG1; 
LULCC = Q_SG2 - Q_SG1; 
CO2 = Q_SG3 - Q_SG2; 

y0 = mean(Clim(year>=1931 & year<=1960,:), 2);
y1 = mean(Clim(year>=1981 & year<=2010,:), 2);
x0 = mean(y0); s0 = std(y0); n0 = length(y0); 
x1 = mean(y1); s1 = std(y1); n1 = length(y1); 
sp = sqrt( ((n0-1)*s0^2 + (n1-1)*s1^2) / (n0+n1-2) );
se = sp * sqrt(1/n0 + 1/n1);
dQ.Var(5) = 'dQ_climate_mstmip';
dQ.dQ(5) = [sprintf('%0.1f ',x1-x0),char(177),sprintf(' %0.1f',1.96*se)];
y = mean(Clim, 2);
mdl = fitlm(1931:2010, y(year>=1931 & year<=2010));
Qbeta.Var(5) = 'Q_climate_mstmip';
Qbeta.Beta(5) = [sprintf('%0.1f ',mdl.Coefficients.Estimate(2)),char(177),sprintf(' %0.1f',1.96*mdl.Coefficients.SE(2))];
for i = 1:nk
    y0 = Clim(year>=1931 & year<=1960,i);
    y1 = Clim(year>=1981 & year<=2010,i);
    x0 = mean(y0); s0 = std(y0); n0 = length(y0); 
    x1 = mean(y1); s1 = std(y1); n1 = length(y1); 
    sp = sqrt( ((n0-1)*s0^2 + (n1-1)*s1^2) / (n0+n1-2) );
    se = sp * sqrt(1/n0 + 1/n1);
    dQens.dQ_clim(i) = [sprintf('%0.1f ',x1-x0),char(177),sprintf(' %0.1f',1.96*se)];
    y = Clim(:,i);
    mdl = fitlm(1931:2010, y(year>=1931 & year<=2010));
    Qbetaens.Beta_clim(i) = [sprintf('%0.1f ',mdl.Coefficients.Estimate(2)),char(177),sprintf(' %0.1f',1.96*mdl.Coefficients.SE(2))];

end

y0 = mean(LULCC(year>=1931 & year<=1960,:), 2);
y1 = mean(LULCC(year>=1981 & year<=2010,:), 2);
x0 = mean(y0); s0 = std(y0); n0 = length(y0); 
x1 = mean(y1); s1 = std(y1); n1 = length(y1); 
sp = sqrt( ((n0-1)*s0^2 + (n1-1)*s1^2) / (n0+n1-2) );
se = sp * sqrt(1/n0 + 1/n1);
dQ.Var(6) = 'dQ_lulcc_mstmip';
dQ.dQ(6) = [sprintf('%0.1f ',x1-x0),char(177),sprintf(' %0.1f',1.96*se)];
y = mean(LULCC, 2);
mdl = fitlm(1931:2010, y(year>=1931 & year<=2010));
Qbeta.Var(6) = 'Q_lulcc_mstmip';
Qbeta.Beta(6) = [sprintf('%0.1f ',mdl.Coefficients.Estimate(2)),char(177),sprintf(' %0.1f',1.96*mdl.Coefficients.SE(2))];
for i = 1:nk
    y0 = LULCC(year>=1931 & year<=1960,i);
    y1 = LULCC(year>=1981 & year<=2010,i);
    x0 = mean(y0); s0 = std(y0); n0 = length(y0); 
    x1 = mean(y1); s1 = std(y1); n1 = length(y1); 
    sp = sqrt( ((n0-1)*s0^2 + (n1-1)*s1^2) / (n0+n1-2) );
    se = sp * sqrt(1/n0 + 1/n1);
    dQens.dQ_lulcc(i) = [sprintf('%0.1f ',x1-x0),char(177),sprintf(' %0.1f',1.96*se)];
    y = LULCC(:,i);
    mdl = fitlm(1931:2010, y(year>=1931 & year<=2010));
    Qbetaens.Beta_lulcc(i) = [sprintf('%0.1f ',mdl.Coefficients.Estimate(2)),char(177),sprintf(' %0.1f',1.96*mdl.Coefficients.SE(2))];

end

y0 = mean(CO2(year>=1931 & year<=1960,:), 2);
y1 = mean(CO2(year>=1981 & year<=2010,:), 2);
x0 = mean(y0); s0 = std(y0); n0 = length(y0); 
x1 = mean(y1); s1 = std(y1); n1 = length(y1); 
sp = sqrt( ((n0-1)*s0^2 + (n1-1)*s1^2) / (n0+n1-2) );
se = sp * sqrt(1/n0 + 1/n1);
dQ.Var(7) = 'dQ_co2_mstmip';
dQ.dQ(7) = [sprintf('%0.1f ',x1-x0),char(177),sprintf(' %0.1f',1.96*se)];
y = mean(CO2, 2);
mdl = fitlm(1931:2010, y(year>=1931 & year<=2010));
Qbeta.Var(7) = 'Q_co2_mstmip';
Qbeta.Beta(7) = [sprintf('%0.1f ',mdl.Coefficients.Estimate(2)),char(177),sprintf(' %0.1f',1.96*mdl.Coefficients.SE(2))];
for i = 1:nk
    y0 = CO2(year>=1931 & year<=1960,i);
    y1 = CO2(year>=1981 & year<=2010,i);
    x0 = mean(y0); s0 = std(y0); n0 = length(y0); 
    x1 = mean(y1); s1 = std(y1); n1 = length(y1); 
    sp = sqrt( ((n0-1)*s0^2 + (n1-1)*s1^2) / (n0+n1-2) );
    se = sp * sqrt(1/n0 + 1/n1);
    dQens.dQ_co2(i) = [sprintf('%0.1f ',x1-x0),char(177),sprintf(' %0.1f',1.96*se)];
    y = CO2(:,i);
    mdl = fitlm(1931:2010, y(year>=1931 & year<=2010));
    Qbetaens.Beta_co2(i) = [sprintf('%0.1f ',mdl.Coefficients.Estimate(2)),char(177),sprintf(' %0.1f',1.96*mdl.Coefficients.SE(2))];

end

writetable(dQ, './output/difference-in-means-instrumental.xlsx');
writetable(dQens, './output/difference-in-means-mstmip-models-instrumental.xlsx');
writetable(Qbeta, './output/trends-instrumental.xlsx');
writetable(Qbetaens, './output/trends-mstmip-models-instrumental.xlsx');
