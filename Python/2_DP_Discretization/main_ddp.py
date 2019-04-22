"""
メインファイル:
状態変数と操作変数を離散化して動的計画法(discretized DP)を解く.
"""

# pylint: disable=C0103

import copy
import math
import time
import numpy as np
#import matplotlib as mpl
import matplotlib.pyplot as plt

## カリブレーション
beta = 0.96 # 割引因子
gamma = 1.0 # 相対的危険回避度(異時点間の代替の弾力性の逆数)
alpha = 0.40 # 資本分配率
delta = 1.00 # 固定資本減耗(0.08)

## 離散化用のパラメータ
nk = 101 # グリッドの数
kmax = 0.5 # 資本グリッドの最大値
#kmax = 10.0 # 資本グリッドの最大値(固定資本減耗=0.08の場合に使用)
kmin = 0.05 # 資本グリッドの最小値(0にすると生産が出来なくなる)

## 収束の基準
it = 1 # ループ・カウンター
maxit = 1000 # 繰り返し計算の最大値
tol = 1.0e-005 # 許容誤差(STEP 2)
dif1 = 1.0 # 価値関数の繰り返し誤差
dif2 = 1.0 # 政策関数の繰り返し誤差

start = time.time()

print("-+- Solve a neoclassical growth model -+-")

## STEP 1(a): グリッド生成

kgrid = np.linspace(kmin, kmax, nk)

## STEP 1(b): 価値関数・政策関数の初期値を設定

vfcn = np.zeros([nk])
pfcn = np.zeros([nk])
Tvfcn = np.zeros([nk])
Tpfcn = np.zeros([nk])
vkp = np.zeros([nk, nk])

## STEP 3: 効用関数の組み合わせ

# 効用関数の初期値 (消費が0以下になる組み合わせにはペナルティ)
util = -10000.0*np.ones([nk, nk])

# 消費が正値になる(k,k')の組み合わせについて効用を計算
for i in range(nk):
    #  あらゆる操作変数k'について:
    for j in range(nk):
        wealth = kgrid[i]**alpha + (1.0-delta)*kgrid[i]
        cons = wealth - kgrid[j]
        if cons > 0:
            if gamma == 1.0:
                util[j, i] = math.log(cons)
            else:
                util[j, i] = cons**(1.0-gamma) / (1.0-gamma)

## STEP 4: 価値関数を繰り返し計算
for it in range(maxit):

    #if dif1 < tol or dif2 < tol:
    #    break
    if dif1 < tol:
        break

    # ベルマン方程式: V(k;k')
    for i in range(nk):
        vkp[0:nk, i] = util[0:nk, i] + beta*vfcn

    # 最適化: 各kについてV(k;k')を最大にするk'を探す
    for i in range(nk):
        Tvfcn[i] = max(vkp[:, i])
        Tpfcn[i] = kgrid[np.argmax(vkp[:, i])]

    # 繰り返し計算誤差を確認
    dif1 = max(abs((Tvfcn-vfcn)/Tvfcn))
    dif2 = max(abs((Tpfcn-pfcn)/Tpfcn))

    # 価値関数・政策関数をアップデート
    # in Python/Julia, "=" does not really assign the "number"
    vfcn = copy.deepcopy(Tvfcn)
    pfcn = copy.deepcopy(Tpfcn)
    print([it, dif1, dif2])

elapsed_time = time.time() - start

print('-+- Computation time -+-')
print(elapsed_time)

## 計算結果をコマンドウィンドウに表示
print('-+- Parameter Values -+-')
print('alpha  beta  gamma      ')
print([alpha, beta, gamma])
print('kmin  kmax  #grid      ')
print([kmin, kmax, nk])

## 解析的解
v_true = np.zeros([nk])
p_true = np.zeros([nk])
AA = (1.0-beta)**(-1) * (math.log(1.0-alpha*beta) + ((alpha*beta)/(1.0-alpha*beta))*math.log(alpha*beta))
BB = alpha/(1.0-alpha*beta)
for i in range(nk):
    v_true[i] = AA + BB*math.log(kgrid[i])
    p_true[i] = beta*alpha*(kgrid[i]**alpha)

## 図を描く
plt.figure()
plt.plot(kgrid, vfcn, linestyle='solid', color='blue', label='approx')
plt.plot(kgrid, v_true, linestyle='dashed', color='red', label='analytical')
plt.title("Value Function")
plt.xlabel("Capital: k")
plt.ylabel("Value Function: V(k)")
plt.legend(loc='lower right')
plt.grid(True)
plt.savefig('Fig3_dndp1.pdf')
plt.show()

plt.figure()
plt.plot(kgrid, pfcn, linestyle='solid', color='blue', label='approx')
plt.plot(kgrid, p_true, linestyle='dashed', color='red', label='analytical')
plt.plot(kgrid, kgrid, linestyle='dotted', color='black', label='45-degree')
plt.title("Policy Function")
plt.xlabel("Current Capital: k")
plt.ylabel("Next Period's Capital: V(k)")
plt.legend(loc='lower right')
plt.grid(True)
plt.savefig('Fig3_dndp2.pdf')
plt.show()
