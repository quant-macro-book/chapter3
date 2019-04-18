%% メインファイル:
% 状態変数と操作変数を離散化して動的計画法(discretized DP)を解く.

clear;
close all;
format short;

%% *** カリブレーション ***
beta  = 0.96; % 割引因子
gamma = 1.0;  % 相対的危険回避度(異時点間の代替の弾力性の逆数)
alpha = 0.40; % 資本分配率
delta = 1.00; % 固定資本減耗(0.08)

% *** 離散化用のパラメータ ***
nk   = 10001;  % グリッドの数
kmax = 0.5;   % 資本グリッドの最大値
%kmax = 10.0; % 資本グリッドの最大値(固定資本減耗=0.08の場合に使用)
kmin = 0.05;  % 資本グリッドの最小値 (0にすると生産が出来なくなる)
%========================

% *** 収束の基準 ***
it = 1;          % ループ・カウンター
maxit = 1000;    % 繰り返し計算の最大値
tol  = 1.0e-005; % 許容誤差(STEP 2)
dif1 = 1;        % 価値関数の繰り返し誤差
dif2 = 1.0;      % 政策関数の繰り返し誤差
count = 1;
%==================

%% 計算開始

tic

disp('')
disp('-+- Solve a neoclassical growth model -+-');

%% STEP 1(a): グリッド生成

kgrid = linspace(kmin, kmax, nk)';

%% STEP 1(b): 価値関数・政策関数の初期値を設定

vfcn  = zeros(nk, 1);
pfcn  = zeros(nk, 1);
Tvfcn = zeros(nk, 1);
Tpfcn = zeros(nk, 1);
vkp   = zeros(nk, nk);
val_tmp = zeros(nk, 4);

%% STEP 3: 効用関数の組み合わせ

% 効用関数の初期値 (消費が0以下になる組み合わせにはペナルティ)
util = -10000.0*ones(nk, nk);

% 消費が正値になる(k,k')の組み合わせについて効用を計算
for i = 1:nk
    %  あらゆる操作変数k'について:
    for j = 1:nk
        wealth = kgrid(i).^alpha + (1.0-delta).*kgrid(i);
        cons = wealth - kgrid(j);
        if cons > 0
           util(j,i) = CRRA(cons, gamma);
        end
    end
end

%% STEP 4: 価値関数を繰り返し計算

while it < maxit && dif1 > tol

    % ベルマン方程式: V(k;k')
    for i = 1:nk
        vkp(:,i) = util(:,i) + beta.*vfcn;
    end
    
    % 最適化: 各kについてV(k;k')を最大にするk'を探す
    [Tvfcn, ploc] = max(vkp);
    Tvfcn = Tvfcn';
    Tpfcn = kgrid(ploc);
    
    % 繰り返し計算誤差を確認
    dif1 = max(abs((Tvfcn-vfcn)./vfcn));
    dif2 = max(abs((Tpfcn-pfcn)./pfcn));

    % 価値関数・政策関数をアップデート
    vfcn = Tvfcn;
    pfcn = Tpfcn;
    fprintf('iteration index: %i, iteration diff of value: %d, iteration diff of policy: %d \n', it, dif1, dif2);

    % 収束途中の繰り返し計算誤差を保存
    % 途中経過を図示する目的なので、通常は不要(むしろ遅くなるので消すべき)
    v_conv(it) = dif1;
    p_conv(it) = dif2;

    % 同じく価値関数の収束を図示する目的で保存(本来は不要)
    if it==1 || it==3 || it==5
       val_tmp(:, count) = vfcn;
       count = count + 1;
    end
    
    it = it + 1;

end

val_tmp(:, 4) = vfcn;

toc

%% 計算結果をコマンドウィンドウに表示

disp('-+- Parameter values -+-');
disp('');
fprintf('beta=%5.2f, gamma=%5.2f, alpha=%5.2f, delta=%5.2f \n', beta, gamma, alpha, delta);
disp(''); 
fprintf('kmin=%5.2f, kmax=%5.2f, #grid=%i \n', kmin, kmax, nk);
disp('');

%% 解析的解

AA = (1.0-beta).^(-1) * (log(1.0-alpha*beta) + ((alpha*beta)/(1.0-alpha*beta))*log(alpha*beta));
BB = alpha/(1.0-alpha*beta);
v_true = AA + BB*log(kgrid);
p_true = beta*alpha*(kgrid.^alpha);

%% オイラー方程式から誤差を測定(非自動化)

nkk = 21;
kgrid2 = linspace(kmin, kmax, nkk)';
pfcn0 = zeros(nkk, 1);
j = 1;

pfcn0(1, 1) = pfcn(1, 1);
pfcn0(2, 1) = pfcn(501, 1);
pfcn0(3, 1) = pfcn(1001, 1);
pfcn0(4, 1) = pfcn(1501, 1);
pfcn0(5, 1) = pfcn(2001, 1);
pfcn0(6, 1) = pfcn(2501, 1);
pfcn0(7, 1) = pfcn(3001, 1);
pfcn0(8, 1) = pfcn(3501, 1);
pfcn0(9, 1) = pfcn(4001, 1);
pfcn0(10, 1) = pfcn(4501, 1);
pfcn0(11, 1) = pfcn(5001, 1);
pfcn0(12, 1) = pfcn(5501, 1);
pfcn0(13, 1) = pfcn(6001, 1);
pfcn0(14, 1) = pfcn(6501, 1);
pfcn0(15, 1) = pfcn(7001, 1);
pfcn0(16, 1) = pfcn(7501, 1);
pfcn0(17, 1) = pfcn(8001, 1);
pfcn0(18, 1) = pfcn(8501, 1);
pfcn0(19, 1) = pfcn(9001, 1);
pfcn0(20, 1) = pfcn(9501, 1);
pfcn0(21, 1) = pfcn(10001, 1);

cons = kgrid2.^alpha + (1.-delta).*kgrid2 - pfcn0(:,1);
LHS  = mu_CRRA(cons, gamma);
kp   = pfcn0(:,1);
kpp  = interp1(kgrid2, pfcn0(:,1), kp);
cons = kp.^alpha + (1.-delta).*kp - kpp;
rent = alpha.*kp.^(alpha-1.0) - delta;
RHS  = beta.*(1.+rent).*mu_CRRA(cons, gamma);
err  = RHS./LHS-1.0;

csvwrite("err_ddp.csv", err);

%% 図を描く

figure;
plot(kgrid, vfcn, '-', 'color', 'blue', 'linewidth', 3); hold('on');
plot(kgrid, v_true, '--', 'color', 'red', 'linewidth', 3); hold('off');
%title('価値関数', 'fontsize', 16);
xlabel('資本保有量：k', 'Fontsize', 16);
ylabel('価値関数：V(k)', 'Fontsize', 16);
xlim([0,kmax]);
legend('近似解', '解析的解', 'Location', 'SouthEast');
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_dndp1.eps','epsc2');

figure;
plot(kgrid, pfcn, '-', 'color', 'blue', 'linewidth', 3); hold('on');
plot(kgrid, p_true, '--', 'color', 'red', 'linewidth', 3);
plot(kgrid, kgrid, ':', 'color', 'black', 'linewidth', 2); hold('off');
%title('政策関数', 'fontsize', 16);
xlabel('今期の資本保有量：k', 'Fontsize', 16);
ylabel("次期の資本保有量：k'", 'Fontsize', 16);
xlim([0,kmax]);
legend('近似解', '解析的解', '45度線', 'Location', 'NorthWest');
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_dndp2.eps','epsc2');

[r, l] = size(v_conv);
iter = linspace(1, l, l)';

figure;
plot(iter, v_conv, '-', 'color', 'blue', 'linewidth', 3);
%title('価値関数の収束', 'fontsize', 16);
xlabel('計算回数', 'Fontsize', 16);
ylabel('繰り返し計算誤差', 'Fontsize', 16);
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_dndp3.eps','epsc2');

figure;
plot(iter, p_conv, '-', 'color', 'blue', 'linewidth', 3);
%title('政策関数の収束', 'fontsize', 16);
xlabel('計算回数', 'Fontsize', 16);
ylabel('繰り返し計算誤差', 'Fontsize', 16);
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_dndp4.eps','epsc2');

figure;
plot(iter, v_conv, '-', 'color', 'blue', 'linewidth', 2); hold('on');
plot(iter, p_conv, ':', 'color', 'red', 'linewidth', 2); hold('off');
%title('価値関数・政策関数の収束', 'fontsize', 16);
xlabel('計算回数', 'Fontsize', 16);
ylabel('繰り返し計算誤差', 'Fontsize', 16);
ylim([0,0.1]);
legend('価値関数', '政策関数', 'Location', 'NorthEast');
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_dndp5.eps','epsc2');

figure;
plot(kgrid, val_tmp(:, 1), '-', 'color', 'blue', 'linewidth', 3); hold('on');
plot(kgrid, val_tmp(:, 2), '--', 'color', 'red', 'linewidth', 3);
plot(kgrid, val_tmp(:, 3), '-.', 'color', 'cyan', 'linewidth', 3);
plot(kgrid, val_tmp(:, 4), ':', 'color', 'magenta', 'linewidth', 3); hold('off');
title('価値関数の収束', 'fontsize', 16);
xlabel('資本保有量：k', 'Fontsize', 16);
ylabel('価値関数：V(k)', 'Fontsize', 16);
xlim([0,kmax]);
legend('it=1', 'it=3', 'it=5', '収束', 'Location', 'East');
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_dndp6.eps','epsc2');

return

