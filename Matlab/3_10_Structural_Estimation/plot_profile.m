clear;
clear global;
close all;
format short;

cprof = readmatrix("consumption_profile.csv");
aprof = readmatrix("asset_profile.csv");

age = linspace(20,105,86);

figure;
plot(age, cprof, '-', 'linewidth', 3);
xlabel('年齢', 'Fontsize', 16);
ylabel('消費', 'Fontsize', 16);
xlim([20,105]);
ylim([0,1.5]);
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_cons_profile.eps','epsc2');

figure;
plot(age, aprof, '-', 'linewidth', 3);
xlabel('年齢', 'Fontsize', 16);
ylabel('資産', 'Fontsize', 16);
xlim([20,105]);
ylim([0,15]);
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_asset_profile.eps','epsc2');

return
