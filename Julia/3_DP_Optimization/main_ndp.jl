#%% 3.4節：Optimizationで無限期間モデルを解く

using Dierckx
using Optim
using Plots
pyplot()

include("GenerateGrid.jl")
include("MyEconFcn.jl")
include("define_struct.jl")
include("calibration.jl")
include("BellmanEq.jl")
include("ndp.jl")

params = calibration()

# 実行
@time vfcn0, pfcn0, dif = ndp(params)

# 最終的な政策関数が得られてから消費関数を計算
wealth = params.kgrid.^params.α + (1-params.δ)*params.kgrid
cfcn = wealth - pfcn0

# 解析解
# AA = (1.0-m.β)^(-1) * (log(1.0-m.α*m.β) + ((m.α*m.β)/(1.0-m.α*m.β))*log(m.α*m.β))
# BB = m.α/(1.0-m.α*m.β)
# v_true = AA .+ BB*log.(m.kgrid)
# p_true = m.α*m.β*(m.kgrid.^m.α);

# オイラー方程式から誤差を測定
# LHS = Utils.mu_CRRA.(cfcn, m.γ)
# kp = pfcn0
# spl = Spline1D(m.kgrid, pfcn0, k = 1, bc = "extrapolate")
# kpp = spl(kp)
# cons = kp.^m.α + (1 - m.δ) * kp - kpp
# rent = m.α*kp.^(m.α - 1) .- m.δ
# RHS = m.β*(1 .+ rent) .* Utils.mu_CRRA.(cons, m.γ)
# err = RHS./LHS .- 1.0;

# プロット

plt = plot(params.kgrid, vfcn0,
    color = :blue,
    legend = :none,
    xlabel = ("現在の資本：k"),
    ylabel = ("価値関数：V(k)"),
    linewidth = 4,
    titlefont = font("HackGen35Nerd", 12),
    guidefont = font("HackGen35Nerd", 12),
    tickfont = font("HackGen35Nerd", 8),
    legend_font_family = ("HackGen35Nerd"),
    legendfontsize = 12,
    framestyle = :semi
)


plt = plot(params.kgrid, pfcn0,
    color = :blue,
    legend = :none,
    xlabel = ("現在の資本：k"),
    ylabel = ("価値関数：V(k)"),
    linewidth = 4,
    titlefont = font("HackGen35Nerd", 12),
    guidefont = font("HackGen35Nerd", 12),
    tickfont = font("HackGen35Nerd", 8),
    legend_font_family = ("HackGen35Nerd"),
    legendfontsize = 12,
    framestyle = :semi
)
