clear;
close all;

addpath('../../Common/')
rng(1000);


cols = GetDefaultColors();

% Load endogenously labelled NPM1 data from Klosin et al. 2020;
[num, txt, raw] = xlsread('Klosin2020_mitosis_NPM1_Fig3F.xlsx');

CC = txt(2:end, 3);
mitosisIdx = strcmp(CC, 'mitosis');

Fl_Mitosis = num(mitosisIdx);
Fl_Interphase = num(~mitosisIdx);

ratioM = mean(Fl_Mitosis) / 7.7;
Fl_Interphase = Fl_Interphase / ratioM;
Fl_Mitosis = Fl_Mitosis / ratioM;

% Load transiently expressed NPM1 data from Riback et al. 2020;
[num, txt, raw] = xlsread('41586_2020_2256_MOESM4_ESM.xlsx', 'Panel b');

Fl_tot = num(:, 1);
Fl_dil = num(:, 2);

p = plot(Fl_tot, Fl_dil, 's', 'Color', 0.6*[1, 1, 1], 'MarkerFaceColor', 0.8*[1, 1, 1]); hold on;
set(p, 'MarkerSize', 6);

p = plot(Fl_Mitosis, Fl_Interphase, '.', 'Color', cols(2, :)); hold on;
set(p, 'MarkerSize', 10);
xlim([0, 100]);
ylim([0, 50]);
yticks([0, 25, 50]);

xlabel('C_{tot} in \muM');
ylabel('C_{dil} in \muM');
legend('Riback et al. 2020', 'Klosin et al. 2020');
col = get(p, 'Color');
box off;

figure;
p = plot(Fl_tot, Fl_dil, 's'); hold on;
set(p, 'MarkerSize', 6, 'Color', 0.6*[1, 1, 1], 'MarkerFaceColor', 0.8*[1, 1, 1]);
p = plot(Fl_Mitosis, Fl_Interphase, '.'); hold on;
set(p, 'MarkerSize', 10, 'Color', cols(2, :));
xlim([0, 15]);
ylim([0, 15]);
plot([0, 15], [0, 15], '--k', 'LineWidth', 2);
xlabel('C_{tot} in \muM');
ylabel('C_{dil} in \muM');
box off;

currAx = gca;
currAx.FontSize = 12;
currAx.XColor = 'k';
currAx.YColor = 'k';

N = length(Fl_Interphase);

idx = randi(N, N, 1000);

Gamma = std(Fl_Mitosis(idx))./mean(Fl_Mitosis(idx)) ./ (std(Fl_Interphase(idx) ./ mean(Fl_Interphase)));

Gamma_Mean = mean(Gamma);
Gamma_Std = std(Gamma);

