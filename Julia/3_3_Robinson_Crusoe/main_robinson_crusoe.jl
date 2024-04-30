#%% 3.3節：ロビンソン・クルーソーモデル

# 必要なパッケージを最初に読み込む
using Plots
pyplot()
using Optim
using Dierckx

include("GenerateGrid.jl")
include("MyEconFcn.jl")
include("define_struct.jl")
include("calibration.jl")
include("BellmanEq.jl")
include("Robinson_crusoe.jl")

# カリブレーションパラメータと変数を読み込む
params = calibration()

# 計算時間を計測
@time vfcn, pfcn, cfcn = Robinson_crusoe(params)

# 解析的解
p_true = zeros(params.nk, params.T)

for t = 1:params.T, i = 1:params.nk
    p_true[i, t] = params.α * params.β * ((1.0-(params.α*params.β)^(params.T-t)) / (1.0-(params.α*params.β)^(params.T-t+1))) * (params.kgrid[i]^params.α)
end

# プロット
plt = plot(params.kgrid, vfcn[:, end],
    legend = :right,
    title = ("各期の価値関数"),
    xlims = (0, 1.1),
    label = ("t=10"),
    xlabel = ("現在の資本：k"),
    ylabel = ("価値関数：V(k)"),
    linewidth = 4,
    titlefont = font("HackGen35Nerd", 12),
    guidefont = font("HackGen35Nerd", 12),
    tickfont = font("HackGen35Nerd", 8),
    framestyle = :semi
)
plot!(params.kgrid, vfcn[:, 9], linewidth =4, label="t=9")
plot!(params.kgrid, vfcn[:, 8], linewidth =4, label="t=8")
plot!(params.kgrid, vfcn[:, 7], linewidth =4, label="t=7")
plot!(params.kgrid, vfcn[:, 1], linewidth =4, label="t=1")
savefig("Fig3_rc1.pdf")


plt = plot(params.kgrid, pfcn[:, end],
    legend = :right,
    title = ("各期の政策関数"),
    xlims = (0, 1.1),
    label = ("t=10"),
    xlabel = ("現在の資本：k"),
    ylabel = ("次期の資本：k'=g(k)"),
    linewidth = 4,
    titlefont = font("HackGen35Nerd", 12),
    guidefont = font("HackGen35Nerd", 12),
    tickfont = font("HackGen35Nerd", 8),
    framestyle = :semi
)
plot!(params.kgrid, pfcn[:, 9], linewidth =4, label="t=9")
plot!(params.kgrid, pfcn[:, 8], linewidth =4, label="t=8")
plot!(params.kgrid, pfcn[:, 7], linewidth =4, label="t=7")
plot!(params.kgrid, pfcn[:, 1], linewidth =4, label="t=1")
savefig("Fig3_rc2.pdf")


plt = plot(params.kgrid, cfcn[:, end],
    legend = :bottomright,
    title = ("各期の消費関数"),
    xlims = (0, 1.1),
    label = ("t=10"),
    xlabel = ("現在の資本：k"),
    ylabel = ("現在の消費：c"),
    linewidth = 4,
    titlefont = font("HackGen35Nerd", 12),
    guidefont = font("HackGen35Nerd", 12),
    tickfont = font("HackGen35Nerd", 8),
    framestyle = :semi
)
plot!(params.kgrid, cfcn[:, 9], linewidth =4, label="t=9")
plot!(params.kgrid, cfcn[:, 8], linewidth =4, label="t=8")
plot!(params.kgrid, cfcn[:, 7], linewidth =4, label="t=7")
plot!(params.kgrid, cfcn[:, 1], linewidth =4, label="t=1")
savefig("Fig3_rc3.pdf")


plt = plot(params.kgrid, pfcn[:, 9],
    legend = :bottomright,
    title = ("真の政策関数と近似した政策"),
    xlims = (0, 1.1),
    label = ("approximation: t=1"),
    xlabel = ("現在の資本：V(k)"),
    ylabel = ("価値関数：k"),
    linewidth = 4,
    titlefont = font("HackGen35Nerd", 12),
    guidefont = font("HackGen35Nerd", 12),
    tickfont = font("HackGen35Nerd", 8),
    framestyle = :semi
)
plot!(params.kgrid, p_true[:, 9], linewidth =4, label="true: t=9")
plot!(params.kgrid, pfcn[:, 1], linewidth =4, label="approximation: t=1")
plot!(params.kgrid, p_true[:, 1], linewidth =4, label="true: t=1")
savefig("Fig3_rc4.pdf")
