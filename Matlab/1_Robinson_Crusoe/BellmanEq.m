function  value = BellmanEq(kprime)
% Function BellmanEq
%  [value] = BellmanEq(kprime)
%
% 目的:
% k'を一つ与えたときのベルマン方程式を返す関数.
% main_robinson_crusoe.mから呼び出して使う.
%
% グローバル変数: beta gamma alpha capital val_fcn kgrid

global beta gamma alpha capital val_fcn kgrid

%% ベルマン方程式

wealth = capital.^alpha;

cons = wealth - kprime;

% 消費が負値の場合、ペナルティを与えてその値が選ばれないようにする
if cons > 0.0
    util = CRRA(cons, gamma);
else
    util = -10000.0;
end

% 次期の価値関数を線形補間
%vnext = interp1(kgrid, val_fcn, kprime, 'linear', 'extrap');

% 次期の価値関数をスプライン補間
vnext = interp1(kgrid, val_fcn, kprime, 'spline');

value = util + beta.*vnext;

%% トリック(1): k'は正の値しか取らないので、ペナルティを与えてその値が選ばれないようにする
if kprime < 0
    value = -1000000.0;
end

%% トリック(2): "最小化"をするので符号を反転
value = -1*value;
 
return
