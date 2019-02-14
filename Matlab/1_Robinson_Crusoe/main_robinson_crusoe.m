%% メインファイル:
% 最適化(optimization)と内挿法(interpolation)をつかってロビンソン・クルーソー経済を解く.

clear;
clear global;
close all;
format short;

global beta gamma alpha capital val_fcn kgrid

%% *** カリブレーション ***
beta  = 0.96; % 割引因子
gamma = 1.0;  % 相対的危険回避度(異時点間の代替の弾力性の逆数)
alpha = 0.4;  % 資本分配率

% *** 離散化用のパラメータ ***
nk   = 11;    % グリッドの数
kmax =  1.0;  % 資本グリッドの最大値
kmin =  0.05; % 資本グリッドの最小値
%==========================

% *** 変数を定義 ***
TT   = 10;            % 無人島に滞在する期間
vfcn = zeros(nk, TT); % 価値関数
pfcn = zeros(nk, TT); % 政策関数
cfcn = zeros(nk, TT); % 消費関数
%=================

%% 計算開始

tic

disp(' ');
fprintf('Solve Robinson Crusoe economy with %d periods \n', TT);

%% グリッドポイントを計算

kgrid = linspace(kmin, kmax, nk)';

%% 最終期(全てを消費)

pfcn(:, TT) = 0; % 全て消費するので貯蓄はゼロ
cfcn(:, TT) = kgrid.^alpha; % 生産=消費
vfcn(:, TT) = CRRA(cfcn(:, TT), gamma); % 消費から得られる効用

%% メインループ

for t = TT-1:-1:1

    fprintf('period %d: \n', t);

    for i = 1:nk

        % グローバル変数を設定
        % fminsearchで使う関数(BellmanEq)に最適化する変数"以外"の変数を渡す(グローバル変数を使わない方法もあるはず)
        capital = kgrid(i);
        val_fcn = vfcn(:, t+1);

        % MATLABの最適化関数(fminsearch)を使ってグリッド上で価値関数と政策関数の値を探す
        % 初期値は0.01
        [pfcn(i, t), vfcn(i, t)] = fminsearch(@BellmanEq, 0.01);

    end

    % 消費関数を計算(必ずしも計算する必要はない)
    cfcn(:, t) = kgrid.^alpha - pfcn(:, t);

    % fminsearchは最小値を探す関数なので符号を反転させる
    vfcn(:, t) = -1*vfcn(:, t);

end

%% 計算結果をコマンドウィンドウに表示

disp('-+- PARAMETER VALUES -+-');
disp('');
fprintf('beta=%5.2f, gamma=%5.2f, alpha=%5.2f \n', beta, gamma, alpha);
disp(''); 
fprintf('kmin=%5.2f, kmax=%5.2f, #grid=%5.2f \n', kmin, kmax, nk);
disp('');

toc

%% 解析的解

p_true = zeros(nk, TT);

for t = 1:TT
    for i = 1:nk
        p_true(i, t) = alpha*beta*( (1-(alpha*beta)^(TT-t)) / (1-(alpha*beta)^(TT-t+1)) )*(kgrid(i)^alpha);
    end
end

%% 図を描く

figure;
plot(kgrid, vfcn(:,TT), '-', 'linewidth', 3); hold('on');
plot(kgrid, vfcn(:,9), '-.', 'linewidth', 3);
plot(kgrid, vfcn(:,8), '--', 'linewidth', 3);
plot(kgrid, vfcn(:,7), ':', 'linewidth', 3);
plot(kgrid, vfcn(:,1), '-', 'linewidth', 3); hold('off');
%title('価値関数', 'fontsize', 16);
xlabel('資本保有量：k_{t}', 'Fontsize', 16);
ylabel('価値関数：V_{t}(k_{t})', 'Fontsize', 16);
legend('10期', '9期', '8期', '7期', '1期', 'Location', 'SouthEast');
grid on;
set(gca,'Fontsize',16);
saveas (gcf,'Fig3_rc1.eps','epsc2');

figure;
plot(kgrid, pfcn(:,TT), '-', 'linewidth', 3); hold('on');
plot(kgrid, pfcn(:,9), '-.', 'linewidth', 3);
plot(kgrid, pfcn(:,8), '--', 'linewidth', 3);
plot(kgrid, pfcn(:,7), ':', 'linewidth', 3);
plot(kgrid, pfcn(:,1), '-', 'linewidth', 3); hold('off');
%title('政策関数', 'fontsize', 16);
xlabel('t期の資本保有量：k_{t}', 'Fontsize', 16);
ylabel("t+1期の資本保有量：k_{t+1}", 'Fontsize', 16);
legend('10期', '9期', '8期', '7期', '1期', 'Location', 'SouthEast');
grid on;
set(gca,'Fontsize',16);
saveas (gcf,'Fig3_rc2.eps','epsc2');

figure;
plot(kgrid, cfcn(:,TT), '-', 'linewidth', 3); hold('on');
plot(kgrid, cfcn(:,9), '-.', 'linewidth', 3);
plot(kgrid, cfcn(:,8), '--', 'linewidth', 3);
plot(kgrid, cfcn(:,7), ':', 'linewidth', 3);
plot(kgrid, cfcn(:,1), '-', 'linewidth', 3); hold('off');
%title('消費関数', 'fontsize', 16);
xlabel('資本保有量：k_{t}', 'Fontsize', 16);
ylabel('消費：c_{t}', 'Fontsize', 16);
legend('10期', '9期', '8期', '7期', '1期', 'Location', 'SouthEast');
grid on;
set(gca,'Fontsize',16);
saveas (gcf,'Fig3_rc3.eps','epsc2');

figure;
plot(kgrid, pfcn(:,9), '-', 'linewidth', 3); hold('on');
plot(kgrid, p_true(:,9), '-.', 'linewidth', 3);
plot(kgrid, pfcn(:,1), ':', 'linewidth', 3);
plot(kgrid, p_true(:,1), '-', 'linewidth', 3); hold('off');
%title('真の政策関数と近似解', 'fontsize', 16);
xlabel('t期の資本保有量：k_{t}', 'Fontsize', 16);
ylabel('t+1期の資本保有量：k_{t+1}', 'Fontsize', 16);
legend('9期の近似解', '9期の解析的解', '1期の近似解', '1期の解析的解', 'Location', 'SouthEast');
grid on;
set(gca,'Fontsize',16);
saveas (gcf,'Fig3_rc4.eps','epsc2');

return
