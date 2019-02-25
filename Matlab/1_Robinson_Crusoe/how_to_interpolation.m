%% main_robinson_crusoe.mを走らせた後で使う:
% 価値関数の内挿をプロットする.
% main_robinson_crusoe.mで使った変数を使うので、先にメインファイルを回す必要あり.

close all;

%% DATA POINTS

figure;
plot(kgrid, vfcn(:, TT-1), 'o', 'color', 'blue', 'MarkerSize', 12, 'linewidth', 3);
xlabel('資本保有量：k', 'Fontsize', 16);
ylabel('グリッド上の価値関数：V_{T-1}(k^{i})', 'Fontsize', 16);
xlim([0, 1.5]);
ylim([-3, 0]);
grid on;
set(gca,'Fontsize',16);
saveas (gcf,'Fig3_data.eps','epsc2');

%% MATLAB関数を使った内挿

nkk = 1001;
kkmin = 0.0;
kkmax = 1.5;
kkgrid = linspace(kkmin, kkmax, nkk);

v_ln = zeros(nkk, 1);
v_cs = zeros(nkk, 1);

for i = 1:nkk
    v_ln(i, 1) = interp1(kgrid, vfcn(:, TT-1), kkgrid(i), 'linear', 'extrap');
    v_cs(i, 1) = interp1(kgrid, vfcn(:, TT-1), kkgrid(i), 'spline');
end

figure;
plot(kkgrid, v_ln, '-', 'color', 'blue', 'linewidth', 3); hold('on');
plot(kkgrid, v_cs, '--', 'color', 'red', 'linewidth', 3); hold('off');
xlabel('資本保有量：k', 'Fontsize', 16);
ylabel('価値関数：V_{T-1}(k)', 'Fontsize', 16);
legend('線形近似', 'スプライン近似', 'Location', 'SouthEast');
xlim([0, 1.5]);
ylim([-3, 0]);
grid on;
set(gca,'Fontsize',16);
saveas (gcf,'Fig3_interp.eps','epsc2');

return;
