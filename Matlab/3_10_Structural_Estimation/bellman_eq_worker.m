function  value = bellman_eq_worker(aprime)
% Function bellman_eq_worker
%  [value] = bellman_eq_worker(kprime)
%
% 目的:
% a'を一つ与えたときのベルマン方程式を返す関数.
% main_lifecycle.mから呼び出して使う.
%
% グローバル変数: beta gamma alpha delta A tran capital vfcn kgrid

global beta gamma asset ss r tran eta endow surv ret_age age zt na nz vfcn_yng vfcn_old agrid

%% ベルマン方程式

% 予算制約と消費
wealth = (1+r)*asset + eta(age)*endow(zt);
cons = wealth - aprime;

% 消費が負値の場合、ペナルティを与えてその値が選ばれないようにする
if cons > 0.0
    util = CRRA(cons, gamma);
else
    util = -10000.0;
end

% 次期の期待効用
if age == ret_age
    % 引退直前であれば次期にはzの不確実性は存在しないので期待値計算は不要
    vnext = interp1(agrid, vfcn_old(:, 1), aprime, 'spline');
    vnext = surv(age)*beta*vnext;
else
    vnext_z = zeros(nz, 1);
    % それぞれの次期のzに対応したvalueを計算
    for zz = 1:nz
        if aprime <= agrid(na)
            vnext_z(zz) = interp1(agrid, vfcn_yng(:, zz, age+1), aprime, 'spline');
        else
            vnext_z(zz) = interp1(agrid, vfcn_yng(:, zz, age+1), aprime, 'linear');
        end
    end
    % 遷移確率行列で期待値を計算
    vnext = surv(age)*beta*dot(tran(zt, :), vnext_z);
end

% 価値関数
value = util + vnext;

%% トリック(1): k'は正の値しか取らないので、ペナルティを与えてその値が選ばれないようにする
if aprime < 0
    value = -1000000.0;
end

%% トリック(2): "最小化"をするので符号を反転
value = -1.0 * value;
 
return
