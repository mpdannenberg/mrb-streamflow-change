% TWO ISSUES HERE (and fix in test_variance... file as well)
    % 1) square R before averaging? revisit...
    % 2) multiple by number of models to get Ns in MsTMIP section


dQ = table('Size',[7 2],'VariableTypes',{'string','string'},...
    'VariableNames',{'Var','dQ'});
dQens = table('Size',[10 4],'VariableTypes',{'string','string','string','string'},...
    'VariableNames',{'Model','dQ_clim','dQ_lulcc','dQ_co2'});

% gage - no 1901-1930 gage data
dQ.Var(1) = 'dQ_gage';
dQ.dQ(1) = 'NA';

% WBM - no 1901-1905 data (burn-in period)
dQ.Var(2) = 'dQ_climate_nat';
dQ.dQ(2) = 'NA';

dQ.Var(3) = 'dQ_climate_anth';
dQ.dQ(3) = 'NA';

dQ.Var(4) = 'dQ_climate_wbm';
dQ.dQ(4) = 'NA';

% MsTMIP
load('./output/MsTMIP-contributions.mat')
dQens.Model = models'; nk = length(models);
Clim = Q_SG1 - Q_RG1; % r_clim = corr(Clim,'rows','pairwise'); r_clim = mean(r_clim(~eye(size(r_clim))));
LULCC = Q_SG2 - Q_SG1; % r_lulcc = corr(LULCC,'rows','pairwise');  r_lulcc = mean(r_lulcc(~eye(size(r_lulcc))));
CO2 = Q_SG3 - Q_SG2; % r_co2 = corr(CO2,'rows','pairwise');  r_co2 = mean(r_co2(~eye(size(r_co2))));

y1 = mean(Clim(year>=1981 & year<=2010,:), 2);
x1 = mean(y1); s1 = std(y1); n1 = length(y1); % * (1 - r_clim^2);
se = s1 / sqrt(n1);
dQ.Var(5) = 'dQ_climate_mstmip';
dQ.dQ(5) = [sprintf('%0.1f ',x1),char(177),sprintf(' %0.1f',1.96*se)];
for i = 1:nk
    y1 = Clim(year>=1981 & year<=2010,i);
    x1 = mean(y1); s1 = std(y1); n1 = length(y1); % * (1 - r_clim^2);
    se = s1 / sqrt(n1);
    dQens.dQ_clim(i) = [sprintf('%0.1f ',x1),char(177),sprintf(' %0.1f',1.96*se)];
end

y1 = mean(LULCC(year>=1981 & year<=2010,:), 2);
x1 = mean(y1); s1 = std(y1); n1 = length(y1); % * (1 - r_lulcc^2);
se = s1 / sqrt(n1);
dQ.Var(6) = 'dQ_lulcc_mstmip';
dQ.dQ(6) = [sprintf('%0.1f ',x1),char(177),sprintf(' %0.1f',1.96*se)];
for i = 1:nk
    y1 = LULCC(year>=1981 & year<=2010,i);
    x1 = mean(y1); s1 = std(y1); n1 = length(y1); % * (1 - r_clim^2);
    se = s1 / sqrt(n1);
    dQens.dQ_lulcc(i) = [sprintf('%0.1f ',x1),char(177),sprintf(' %0.1f',1.96*se)];
end

y1 = mean(CO2(year>=1981 & year<=2010,:), 2);
x1 = mean(y1); s1 = std(y1); n1 = length(y1); % * (1 - r_co2^2);
se = s1 / sqrt(n1);
dQ.Var(7) = 'dQ_co2_mstmip';
dQ.dQ(7) = [sprintf('%0.1f ',x1),char(177),sprintf(' %0.1f',1.96*se)];
for i = 1:nk
    y1 = CO2(year>=1981 & year<=2010,i);
    x1 = mean(y1); s1 = std(y1); n1 = length(y1); % * (1 - r_clim^2);
    se = s1 / sqrt(n1);
    dQens.dQ_co2(i) = [sprintf('%0.1f ',x1),char(177),sprintf(' %0.1f',1.96*se)];
end

writetable(dQ, './output/difference-in-means-preindustrial.xlsx');
writetable(dQens, './output/difference-in-means-mstmip-models-preindustrial.xlsx');
