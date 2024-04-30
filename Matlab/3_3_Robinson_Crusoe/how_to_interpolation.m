%% main_robinson_crusoe.mを走らせた後で使う:
% 価値関数の内挿をプロットする.
% main_robinson_crusoe.mで使った変数を使うので、先にメインファイルを回す必要あり.

close all;


%% Jupyter Notebookの例

x = linspace(-10, 10, 5);
xx = linspace(-15, 15, 1001);

% パラメータの設定
a = 0.75;
b = 2.0;
c = -10.0;

y = test_function(x, a, b, c);

figure;
plot(x, y, 'o', 'MarkerSize', 12, 'linewidth', 3);
xlabel('x', 'Fontsize', 16);
ylabel('y = f(x)', 'Fontsize', 16);
xlim([-15, 15]);
grid on;
set(gca,'Fontsize',16);
saveas (gcf,'Fig3_data2.eps','epsc2');

figure;
plot(x, y, 'o', 'color', 'black', 'MarkerSize', 12, 'linewidth', 3);
xlabel('x', 'Fontsize', 16);
ylabel('y = f(x)', 'Fontsize', 16);
xlim([-15, 15]);
grid on;
set(gca,'Fontsize',16);
saveas (gcf,'Fig3_data2_bk.eps','epsc2');

y_ln = zeros(1, 1001);
y_sp = zeros(1, 1001);

for i = 1:1001
    y_ln(i) = interp1(x, y, xx(i), 'linear', 'extrap');
    y_sp(i) = spline(x, y, xx(i));
end

figure;
plot(xx, y_ln, '-', 'linewidth', 3); hold('on');
plot(xx, y_sp, '--', 'linewidth', 3); hold('off');
xlabel('x', 'Fontsize', 16);
ylabel('y = f(x)', 'Fontsize', 16);
legend('線形近似', 'スプライン近似', 'Location', 'SouthEast');
xlim([-15, 15]);
grid on;
set(gca,'Fontsize',16);
saveas (gcf,'Fig3_interp2.eps','epsc2');

figure;
plot(xx, y_ln, '-', 'color', 'black', 'linewidth', 3); hold('on');
plot(xx, y_sp, '--', 'color', 'black', 'linewidth', 3); hold('off');
xlabel('x', 'Fontsize', 16);
ylabel('y = f(x)', 'Fontsize', 16);
legend('線形近似', 'スプライン近似', 'Location', 'SouthEast');
xlim([-15, 15]);
grid on;
set(gca,'Fontsize',16);
saveas (gcf,'Fig3_interp2_bk.eps','epsc2');


%% DATA POINTS

figure;
plot(kgrid, vfcn(:, TT-1), 'o', 'MarkerSize', 12, 'linewidth', 3);
xlabel('資本保有量：k', 'Fontsize', 16);
ylabel('グリッド上の価値関数：V_{T-1}(k^{i})', 'Fontsize', 16);
xlim([0, 1.5]);
ylim([-3, 0]);
grid on;
set(gca,'Fontsize',16);
saveas (gcf,'Fig3_data.eps','epsc2');

figure;
plot(kgrid, vfcn(:, TT-1), 'o', 'color', 'black', 'MarkerSize', 12, 'linewidth', 3);
xlabel('資本保有量：k', 'Fontsize', 16);
ylabel('グリッド上の価値関数：V_{T-1}(k^{i})', 'Fontsize', 16);
xlim([0, 1.5]);
ylim([-3, 0]);
grid on;
set(gca,'Fontsize',16);
saveas (gcf,'Fig3_data_bk.eps','epsc2');

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
plot(kkgrid, v_ln, '-', 'linewidth', 3); hold('on');
plot(kkgrid, v_cs, '--', 'linewidth', 3); hold('off');
xlabel('資本保有量：k', 'Fontsize', 16);
ylabel('価値関数：V_{T-1}(k)', 'Fontsize', 16);
legend('線形近似', 'スプライン近似', 'Location', 'SouthEast');
xlim([0, 1.5]);
ylim([-3, 0]);
grid on;
set(gca,'Fontsize',16);
saveas (gcf,'Fig3_interp.eps','epsc2');


figure;
plot(kkgrid, v_ln, '-', 'color', 'black', 'linewidth', 3); hold('on');
plot(kkgrid, v_cs, '--', 'color', 'black', 'linewidth', 3); hold('off');
xlabel('資本保有量：k', 'Fontsize', 16);
ylabel('価値関数：V_{T-1}(k)', 'Fontsize', 16);
legend('線形近似', 'スプライン近似', 'Location', 'SouthEast');
xlim([0, 1.5]);
ylim([-3, 0]);
grid on;
set(gca,'Fontsize',16);
saveas (gcf,'Fig3_interp_bk.eps','epsc2');


return;
