%% メインファイル:
% 状態変数のみ離散化して操作変数は連続的に値を取る場合の動的計画法(parametric DP)の解法.
% アルゴリズムの詳細は、Johnson et al. (1993)を参照

clear;
clear global;
close all;
format short;

global beta gamma asset ss r tran eta endow surv ret_age age zt na nz vfcn_yng vfcn_old agrid

%% *** カリブレーション ***
beta  = 0.98; % 割引因子
gamma = 1.0;  % 相対的危険回避度(異時点間の代替の弾力性の逆数)
r     = 0.04; % 利子率
ss    = 0.5;  % 所得代替率でカリブレート

eta = readmatrix("earnings_profiles.csv");
surv = readmatrix("surv.csv");

nz = 3;
endow = [0.8027, 1.0, 1.2457];
tran = [0.7451 0.2528 0.0021; 0.1360 0.7281 0.1360; 0.0021 0.2528 0.7451];

max_age = 86; %20歳から105歳まで
ret_age = 45; %20歳から64歳まで

% *** 離散化用のパラメータ ***
na   = 101;  % グリッドの数
amax = 40.0; % 資産グリッドの最大値
amin = 0.0;  % 資産グリッドの最小値 
%========================

%% 計算開始

tic

disp('')
disp('-+- Solve a life cycle model -+-');

%% STEP 1(a): グリッド生成

agrid = linspace(amin, amax, na)';

%% STEP 1(b): 価値関数・政策関数を入れるための変数

pfcn_yng = zeros(na, nz, ret_age);
vfcn_yng = zeros(na, nz, ret_age);

pfcn_old = zeros(na, max_age - ret_age);
vfcn_old = zeros(na, max_age - ret_age);

%% STEP 4: 価値関数を後ろ向きに計算

% 最終期はすべての資源を使い切る
for i = 1:na
    vfcn_old(i, max_age-ret_age) = CRRA((1+r)*agrid(i) + ss, gamma);
    pfcn_old(i, max_age-ret_age) = 0.0;
end

figure;
plot(agrid, vfcn_old(:, max_age-ret_age), '-', 'linewidth', 3);
xlabel('貯蓄', 'Fontsize', 16);
ylabel('価値関数', 'Fontsize', 16);
grid on;

% 引退後の価値関数
for t = max_age-1:-1:ret_age+1

    fprintf('age index: %i \n', t);

    for i = 1:na
    
        % グローバル変数を設定
        % fminsearchで使う関数(BellmanEq)に最適化する変数"以外"の変数を渡す(グローバル変数を使わない方法もあるはず)
        asset = agrid(i);
        age = t;
    
        % MATLABの最適化関数(fminsearch)を使ってグリッド上で価値関数と政策関数の値を探す
        % 初期値は0.01
        [pfcn_old(i, age-ret_age), vfcn_old(i, age-ret_age)] = fminsearch(@bellman_eq_retired, 0.01);
    
    end

    % fminsearchは最小値を探す関数なので符号を反転させる
    vfcn_old(:, age-ret_age) = -1*vfcn_old(:, age-ret_age);

end

figure;
plot(agrid, vfcn_old(:, 1), '-', 'linewidth', 3);
xlabel('貯蓄', 'Fontsize', 16);
ylabel('価値関数', 'Fontsize', 16);
grid on;

figure;
plot(agrid, pfcn_old(:, 1), '-', 'linewidth', 3);
xlabel('貯蓄', 'Fontsize', 16);
ylabel('政策関数', 'Fontsize', 16);
grid on;

% 労働者の価値関数
for t = ret_age:-1:1

    fprintf('age index: %i \n', t);

    for z = 1:nz
        for i = 1:na
        
            % グローバル変数を設定
            % fminsearchで使う関数(BellmanEq)に最適化する変数"以外"の変数を渡す(グローバル変数を使わない方法もあるはず)
            asset = agrid(i);
            age = t;
            zt  = z;

            % MATLABの最適化関数(fminsearch)を使ってグリッド上で価値関数と政策関数の値を探す
            % 初期値は0.01
            [pfcn_yng(i, z, age), vfcn_yng(i, z, age)] = fminsearch(@bellman_eq_worker, 0.01);
        
        end
    end

    % fminsearchは最小値を探す関数なので符号を反転させる
    vfcn_yng(:, :, age) = -1*vfcn_yng(:, :, age);

end

toc

%% 図を描く

age = linspace(20, 64, 45);

figure;
plot(age, eta, '-', 'linewidth', 3);
xlabel('年齢', 'Fontsize', 16);
ylabel('年齢ごとの労働生産性', 'Fontsize', 16);
xlim([20,65]);
ylim([0,1.5]);
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_eta.eps','epsc2');

age = linspace(20, 105, 86);

figure;
plot(age, surv, '-', 'linewidth', 3);
xlabel('年齢', 'Fontsize', 16);
ylabel('生存確率', 'Fontsize', 16);
xlim([20,105]);
ylim([0,1.0]);
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig3_surv.eps','epsc2');

figure;
plot(agrid, vfcn_yng(:, 1, 1), '-', 'linewidth', 3); hold('on');
plot(agrid, vfcn_yng(:, 2, 1), '-', 'linewidth', 3);
plot(agrid, vfcn_yng(:, 3, 1), '-', 'linewidth', 3); hold('off');
xlabel('貯蓄', 'Fontsize', 16);
ylabel('価値関数', 'Fontsize', 16);
grid on;

figure;
plot(agrid, pfcn_old(:, 1, 1), '-', 'linewidth', 3); hold('on');
plot(agrid, pfcn_old(:, 2, 1), '-', 'linewidth', 3);
plot(agrid, pfcn_old(:, 3, 1), '-', 'linewidth', 3); hold('off');
xlabel('貯蓄', 'Fontsize', 16);
ylabel('政策関数', 'Fontsize', 16);
grid on;

return
