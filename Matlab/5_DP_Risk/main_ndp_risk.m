%% メインファイル:
% 状態変数のみ離散化して操作変数は連続的に値を取る場合の動的計画法(parametric DP)の解法.
% アルゴリズムの詳細は、Johnson et al. (1993)を参照

clear;
clear global;
close all;
format short;

global beta gamma alpha delta A tran capital vfcn kgrid

%% *** カリブレーション ***
beta  = 0.96; % 割引因子
gamma = 1.0;  % 相対的危険回避度(異時点間の代替の弾力性の逆数)
alpha = 0.40; % 資本分配率
delta = 1.00; % 固定資本減耗(0.08)
tfp   = [1.01, 0.99];
Pi_A  = [1-0.125 0.125; 0.125 1-0.125];

% *** 離散化用のパラメータ ***
nk   = 21;    % グリッドの数
kmax = 0.5;   % 資本グリッドの最大値
%kmax = 10.0; % 資本グリッドの最大値(固定資本減耗=0.08の場合に使用)
kmin = 0.05;  % 資本グリッドの最小値 (0にすると生産が出来なくなる)
%========================

% *** 収束の基準 ***
it = 1;          % ループ・カウンター
maxit = 1000;    % 繰り返し計算の最大値
tol  = 1.0e-005; % 許容誤差(STEP 2)
dif1 = 1.0;      % 価値関数の繰り返し誤差
dif2 = 1.0;      % 政策関数の繰り返し誤差
count = 1;
%=================

%% 計算開始

tic

disp('')
disp('-+- Solve a neoclassical growth model with TFP shock -+-');

%% STEP 1(a): グリッド生成

kgrid = linspace(kmin, kmax, nk)';
% kgrid = grid_exp1(kmin, kmax, nk)';
% kgrid = grid_exp2(kmin, kmax, nk)';
% kgrid = grid_exp3(kmin, kmax, nk)';

%% STEP 1(b): 価値関数・政策関数の初期値を当て推量

% Aの状態はgoodとbadの2つ
pfcn0 = zeros(nk, 2);
vfcn0 = zeros(nk, 2);
vfcn0(:, 1) = CRRA(tfp(1)*kgrid.^alpha + (1.-delta).*kgrid, gamma);
vfcn0(:, 2) = CRRA(tfp(2)*kgrid.^alpha + (1.-delta).*kgrid, gamma);

pfcn1 = zeros(nk, 2);
vfcn1 = zeros(nk, 2);

%% STEP 4: 価値関数を繰り返し計算

while it < maxit && dif1(1) > tol

    fprintf('iteration index: %i \n', it);
    fprintf('value function iteration error: %e\n', dif1);
    fprintf('policy function iteration error: %e\n', dif2);

    for z = 1:2
        for i = 1:nk
    
            % グローバル変数を設定
            % fminsearchで使う関数(BellmanEq)に最適化する変数"以外"の変数を渡す(グローバル変数を使わない方法もあるはず)
            capital = kgrid(i);
            A = tfp(z);
            tran = Pi_A(z, :);
            vfcn = vfcn0;
    
            % MATLABの最適化関数(fminsearch)を使ってグリッド上で価値関数と政策関数の値を探す
            % 初期値は0.01
            [pfcn1(i,z), vfcn1(i,z)] = fminsearch(@BellmanEq, 0.01);
    
        end
    end

    % fminsearchは最小値を探す関数なので符号を反転させる
    vfcn1 = -1*vfcn1;

    % 繰り返し計算誤差を確認
    dif1 = max(abs((vfcn1-vfcn0)./vfcn0));
    dif2 = max(abs((pfcn1-pfcn0)./pfcn0));

    % 価値関数・政策関数をアップデート
    vfcn0 = vfcn1;
    pfcn0 = pfcn1;

    it = it + 1;

end

%% 計算結果をコマンドウィンドウに表示

disp('-+- PARAMETER VALUES -+-');
disp('');
fprintf('beta=%5.2f, gamma=%5.2f, alpha=%5.2f, delta=%5.2f \n', beta, gamma, alpha, delta);
disp(''); 
fprintf('kmin=%5.2f, kmax=%5.2f, #grid=%i \n', kmin, kmax, nk);
disp('');

toc

%% 図を描く

figure;
plot(kgrid, vfcn0(:, 1), '-', 'linewidth', 3); hold('on');
plot(kgrid, vfcn0(:, 2), '-.', 'linewidth', 3); hold('off');
%title('価値関数', 'fontsize', 16);
xlabel('資本保有量：k', 'Fontsize', 16);
ylabel('価値関数：V(k)', 'Fontsize', 16);
xlim([0, kmax]);
legend('好況', '不況', 'Location', 'NorthWest');
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_pndp_risk1.eps','epsc2');

figure;
plot(kgrid, pfcn0(:, 1), '-', 'linewidth', 3); hold('on');
plot(kgrid, pfcn0(:, 2), '-.', 'linewidth', 3);
plot(kgrid, kgrid, ':', 'linewidth', 2); hold('off');
%title('政策関数', 'fontsize', 16);
xlabel('今期の資本保有量：k', 'Fontsize', 16);
ylabel("次期の資本保有量：k'", 'Fontsize', 16);
xlim([0, kmax]);
ylim([0, kmax]);
legend('好況', '不況', '45度線', 'Location', 'NorthWest');
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_pndp_risk2.eps','epsc2');


%% 白黒の図

figure;
plot(kgrid, vfcn0(:, 1), '-', 'color', 'black', 'linewidth', 3); hold('on');
plot(kgrid, vfcn0(:, 2), '-.', 'color', 'black', 'linewidth', 3); hold('off');
%title('価値関数', 'fontsize', 16);
xlabel('資本保有量：k', 'Fontsize', 16);
ylabel('価値関数：V(k)', 'Fontsize', 16);
xlim([0, kmax]);
legend('好況', '不況', 'Location', 'NorthWest');
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_pndp_risk1_bk.eps','epsc2');

figure;
plot(kgrid, pfcn0(:, 1), '-', 'color', 'black', 'linewidth', 3); hold('on');
plot(kgrid, pfcn0(:, 2), '-.', 'color', 'black', 'linewidth', 3);
plot(kgrid, kgrid, ':', 'color', 'black', 'linewidth', 2); hold('off');
%title('政策関数', 'fontsize', 16);
xlabel('今期の資本保有量：k', 'Fontsize', 16);
ylabel("次期の資本保有量：k'", 'Fontsize', 16);
xlim([0, kmax]);
ylim([0, kmax]);
legend('好況', '不況', '45度線', 'Location', 'NorthWest');
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_pndp_risk2_bk.eps','epsc2');

return
