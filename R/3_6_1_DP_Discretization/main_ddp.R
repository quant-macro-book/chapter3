# メインファイル
# 状態変数と操作変数を離散化して動的計画法(discretized DP)を解く.

# environmentから変数をクリア
rm(list=ls())

# カリブレーション
beta  <- 0.96 # 割引因子
gamma <- 1.0  # 相対的危険回避度(異時点間の代替の弾力性の逆数)
alpha <- 0.40 # 資本分配率
delta <- 1.00 # 固定資本減耗

# 離散化用のパラメータ
nk   <- 1001  # グリッドの数
kmax <- 0.5   # 資本グリッドの最大値
kmin <- 0.05  # 資本グリッドの最小値 (0にすると生産が出来なくなる)

# 収束の基準
maxit <- 1000    # 繰り返し計算の最大値
tol   <- 1.0e-005 # 許容誤差(STEP 2)
dif1  <- 1.0      # 価値関数の繰り返し誤差
dif2  <- 1.0      # 政策関数の繰り返し誤差
#----------------------------------------

# 計算開始
time_start <- proc.time()

cat('', "\n")
cat('-+- Solving an optimal growth model -+-', "\n")
cat('', "\n")

# STEP 1(a):グリッド生成
kgrid <- seq(from = kmin, to = kmax, length.out = nk)

# STEP 1(b):価値関数・政策関数の初期値を設定
vfcn  <- matrix(data = 0.0, nrow = nk, ncol = 1)
pfcn  <- matrix(data = 0.0, nrow = nk, ncol = 1)
Tvfcn <- matrix(data = 0.0, nrow = nk, ncol = 1)
Tpfcn <- matrix(data = 0.0, nrow = nk, ncol = 1)
vkp   <- matrix(data = 0.0, nrow = nk, ncol = nk)

# STEP 3: 効用関数の組み合わせ

# 効用関数の初期値 (消費が0以下になる組み合わせにはペナルティ)
util <- matrix(data = -10000.0, nrow = nk, ncol = nk)

# 消費が正値になる(k,k')の組み合わせについて効用を計算
for (i in 1:nk) {
    #  あらゆる操作変数k'について:
    for (j in 1:nk) {
        wealth <- kgrid[i]^alpha + (1.0-delta)*kgrid[i]
        cons   <- wealth - kgrid[j]
        if (cons > 0) {
            if (gamma == 1.0) {
                util[j, i] <- log(cons)
            } else {
                util[j, i] <- cons^(1.0-gamma) / (1.0-gamma)
            }
        }
    }
}

# STEP 4: 価値関数を繰り返し計算
for (it in 1:maxit) {

    #if (dif1 < tol || dif2 < tol) break
    if (dif1 < tol) break

    # ベルマン方程式: V(k;k')
    for (i in 1:nk) {
        vkp[1:nk, i] = util[1:nk, i] + beta*vfcn
    }

    # 最適化: 各kについてV(k;k')を最大にするk'を探す
    for (i in 1:nk) {
        Tvfcn[i] <- max(vkp[1:nk, i])
        Tpfcn[i] <- kgrid[which.max(vkp[1:nk, i])]
    }

    # 繰り返し計算誤差を確認
    dif1 <- max(abs((Tvfcn-vfcn)/vfcn))
    dif2 <- max(abs((Tpfcn-pfcn)/pfcn))

    # 価値関数・効用関数をアップデート
    vfcn <- Tvfcn
    pfcn <- Tpfcn
    cat(" iteration index:", it, ", iter dif of value:", dif1, ", iter dif of policy:", dif2, "\n")

}

# 計算時間をカウント終了
cat(" Time = ", proc.time() - time_start, "\n")
cat('', "\n")

# 計算結果を表示
cat('-+- Parameter values -+-')
cat(" beta = ", beta, ", gamma = ", gamma, ", alpha = ", alpha, ", delta = ", delta, "\n")
cat(" kmin:", kmin, ", kmax:", kmax, ", #grid:", nk, "\n")
cat("\n")

# 解析的解
AA = (1.0-beta)^(-1) * (log(1.0-alpha*beta) + ((alpha*beta)/(1.0-alpha*beta))*log(alpha*beta))
BB = alpha/(1.0-alpha*beta)
v_true = AA + BB*log(kgrid)
p_true = beta*alpha*(kgrid**alpha)

# 図を描く

library("ggplot2")
library("extrafont")

# Figure 1: value function
df <- data.frame(k = kgrid, v1 = vfcn, v2 = v_true)
df1 <- data.frame(cap = df$k, v1 = df$v1, group="Approx")
df2 <- data.frame(cap = df$k, v2 = df$v2, group="Analytical")

g1 <- ggplot()
g1 <- g1 + geom_line(data=df1, aes(x = cap, y = v1, color=group), size=1.5)
g1 <- g1 + geom_line(data=df2, aes(x = cap, y = v2, color=group), size=1.5, linetype="dashed")
g1 <- g1 + scale_color_manual(name="Legend", values=c("Approx"="blue", "Analytical"="red"))
g1 <- g1 + labs(title="Value Function", x="Capital: k", y="Value Function: V(k)")
plot(g1)
ggsave(file = "Fig_dndp1.eps", width = 8, height = 6)

# Figure 2: policy function
fp <- data.frame(k = kgrid, v1 = pfcn, v2 = p_true, v3 = kgrid)
fp1 <- data.frame(cap = fp$k, v1 = fp$v1, group="Approx")
fp2 <- data.frame(cap = fp$k, v2 = fp$v2, group="Analytical")
fp3 <- data.frame(cap = fp$k, v3 = fp$v3, group="45-degree")

g2 <- ggplot()
g2 <- g2 + geom_line(data=fp1, aes(x = cap, y = v1, color=group), size=1.5)
g2 <- g2 + geom_line(data=fp2, aes(x = cap, y = v2, color=group), size=1.5)
g2 <- g2 + geom_line(data=fp3, aes(x = cap, y = v3, color=group), size=1.5, linetype="dotted")
g2 <- g2 + scale_color_manual(name="Legend", values=c("Approx"="blue", "Analytical"="red", "45-degree"="black"))
g2 <- g2 + labs(title="Policy Function", x="Current Capital", y="Next Period's Capital")
plot(g2)
ggsave(file = "Fig_dndp2.eps", width = 8, height = 6)
