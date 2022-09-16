#%% 3.4節：グリッドサーチで無限期間モデルを解く

using Plots
pyplot()

include("GenerateGrid.jl")
include("MyEconFcn.jl")
include("define_struct.jl")
include("calibration.jl")
include("ddp.jl")

params = calibration()

# 実行
@time vfcn, pfcn, dif, val_tmp = ddp(params)

# 政策関数
policy = zeros(params.nk)
for i = 1:params.nk
    policy[i] = params.kgrid[pfcn[i]]
end

# 収束した価値関数
util = zeros(params.nk)
valfn = zeros(params.nk)

for i = 1:params.nk
    cons = params.kgrid[i]^params.α + (1-params.δ)*params.kgrid[i] - policy[i]
    util[i] = MyEconFcn.crra(cons, params.γ)
    valfn[i] = util[i]/(1-params.β)
end

# 解析解：δ=1の場合のみ
# AA = (1.0-m.β)^(-1) * (log(1.0-m.α*m.β) + ((m.α*m.β)/(1.0-m.α*m.β))*log(m.α*m.β))
# BB = m.α/(1.0-m.α*m.β)
# v_true = AA .+ BB*log.(m.kgrid)
# p_true = m.α*m.β*(m.kgrid.^m.α);


# オイラー方程式から誤差を測定
# nkk = 21
# kgrid2 = collect(LinRange(m.kgrid[1], m.kgrid[end], nkk))
# pfcn0 = zeros(nkk)

# for i in 1:nkk
#     pfcn0[i] = pfcn[500*(i-1)+1]
# end

# cons = kgrid2.^m.α + (1.0 - m.δ)*kgrid2 - pfcn0

# LHS = Utils.mu_CRRA.(cons, m.γ)
# kp = pfcn0
# spl = Spline1D(kgrid2, pfcn0, k = 1, bc = "extrapolate")
# kpp = spl(kp)
# cons = kp.^m.α + (1 - m.δ) * kp - kpp
# rent = m.α*kp.^(m.α - 1) .- m.δ
# RHS = m.β*(1 .+ rent) .* Utils.mu_CRRA.(cons, m.γ)
# err = RHS./LHS .- 1.0;

# プロット
plt = plot(params.kgrid, valfn,
    color = :blue,
    legend = :bottomright,
    label = ("収束"),
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


plt = plot(params.kgrid, val_tmp[:, 1],
    color = :blue,
    legend = :bottomright,
    label = ("t=1"),
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
plot!(params.kgrid, val_tmp[:, 2], linewidth =4, label="t=3")
plot!(params.kgrid, val_tmp[:, 3], linewidth =4, label="t=5")


plt = plot(params.kgrid, policy,
    color = :blue,
    legend = :none,
    label = ("収束"),
    xlabel = ("現在の資本：k"),
    ylabel = ("次期の資本：k'=g(k)"),
    linewidth = 4,
    titlefont = font("HackGen35Nerd", 12),
    guidefont = font("HackGen35Nerd", 12),
    tickfont = font("HackGen35Nerd", 8),
    legend_font_family = ("HackGen35Nerd"),
    legendfontsize = 12,
    framestyle = :semi
)
