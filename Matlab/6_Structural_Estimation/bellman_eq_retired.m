function  value = bellman_eq_retired(aprime)
% Function bellman_eq_retired
%  [value] = bellman_eq_retired(kprime)
%
% 目的:
% a'を一つ与えたときのベルマン方程式を返す関数.
% main_lifecycle.mから呼び出して使う.
%
% グローバル変数: beta gamma alpha delta A tran capital vfcn kgrid

global beta gamma asset ss r surv ret_age age vfcn_old agrid

%% ベルマン方程式

% 予算制約と消費
wealth = (1+r)*asset + ss;
cons = wealth - aprime;

% 消費が負値の場合、ペナルティを与えてその値が選ばれないようにする
if cons > 0.0
    util = CRRA(cons, gamma);
else
    util = -10000.0;
end

% 次期の価値関数をスプライン補間
vnext = interp1(agrid, vfcn_old(:, age-ret_age+1), aprime, 'spline');

% 価値関数
value = util + surv(age)*beta*vnext;

%% トリック(1): k'は正の値しか取らないので、ペナルティを与えてその値が選ばれないようにする
if aprime < 0
    value = -1000000.0;
end

%% トリック(2): "最小化"をするので符号を反転
value = -1.0 * value;
 
return
