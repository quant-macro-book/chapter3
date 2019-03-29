using LinearAlgebra
using Dierckx
using Optim
using Plots
using DelimitedFiles
using LaTeXStrings

include("BellmanEq.jl")
include("Models.jl")
include("Utils.jl")
include("ndp.jl")

m = Models(nk=21, kmax= 0.5) # パラメータを含む構造体のインスタンス

# 実行
vfcn0, pfcn0, dif = ndp(m)
@time ndp(m)
#最終的な政策関数が得られてから消費関数を計算
wealth = m.kgrid.^m.α + (1 -m.δ)*m.kgrid
cfcn = wealth - pfcn0;

#解析解
AA = (1.0-m.β)^(-1) * (log(1.0-m.α*m.β) + ((m.α*m.β)/(1.0-m.α*m.β))*log(m.α*m.β))
BB = m.α/(1.0-m.α*m.β)
v_true = AA .+ BB*log.(m.kgrid)
p_true = m.α*m.β*(m.kgrid.^m.α);

# オイラー方程式から誤差を測定
LHS = Utils.mu_CRRA.(cfcn, m.γ)
kp = pfcn0
spl = Spline1D(m.kgrid, pfcn0, k = 1, bc = "extrapolate")
kpp = spl(kp)
cons = kp.^m.α + (1 - m.δ) * kp - kpp
rent = m.α*kp.^(m.α - 1) .- m.δ
RHS = m.β*(1 .+ rent) .* Utils.mu_CRRA.(cons, m.γ)
err = RHS./LHS .- 1.0;

#=プロット=#
plt1 = plot(m.kgrid, vfcn0,
            label ="Approximated solution",
            line = (1.5, :solid), color ="blue",
            legend=:bottomright)
plot!(m.kgrid, v_true,
      label = "True solution",
      line = (1.5, :dash), color ="red")

xlabel!(L"\mathrm{Amount\ of\ capital :  k }")
ylabel!(L"\mathrm{Value\  function  :  V(k)} ")
xlims!(0, m.kgrid[end])
#savefig("Fig3_pndp1.pdf")
display(plt1)

plt2 =plot(m.kgrid, pfcn0,
           label ="Approximated solution",
           line = (1.5, :solid), color ="blue",
           legend=:bottomright)
plot!(m.kgrid, p_true,
      label = "True solution",
      line = (1.5, :dash), color ="red")

plot!(m.kgrid, m.kgrid,
      label = "45 degree line",
      line = (1.5, :dot), color="black")
xlabel!(L"\mathrm{Amount\ of\ capital\ in\ the\ current\ period :  k }")
ylabel!(L"\mathrm{Amount\ of\ capital\ in\ the\ next\ period :  k'}")
xlims!(0, m.kgrid[end])
#savefig("Fig3_pndp2.pdf")
display(plt2)

plt3 = plot(m.kgrid, cfcn, label="", color="blue", line = 1.5)
xlabel!(L"\mathrm{Amount\ of\ capital :  k }")
ylabel!(L"\mathrm{Consumption :  c }")
#savefig("Fig3_pndp3.pdf")
display(plt3)

styles = filter((s->begin
                s in Plots.supported_styles()
            end), [:solid, :dash])
styles = reshape(styles, 1, length(styles))

x = 1:m.maxiter
plt4 = plot(x, dif', line = (1.5, styles), color =["blue" "red"],
     label =["Value function" "Policy function"])
xlims!(0,250)
ylims!(0, 0.1)
ylabel!("Approximation error")
xlabel!("Number of iterations")
#savefig("Fig3_pndp4.pdf")
display(plt4)

plt5 = plot(m.kgrid, err,label ="", color = "blue")
ylims!(-15e-4, 5e-4)
ylabel!("Euler-equation errors")
xlabel!(L"\mathrm{Amount\ of\ capital :  k }")
#savefig("Fig3_pndp5.pdf")
display(plt5)

dir = pwd()
text_path = dir * "/Julia/3_DP_Optimization"
cd(text_path)
err2 = readdlm("err_ddp.txt")
err2 = vec(err2);

plt6 = plot(m.kgrid, err, color = "red",line =(1.5, :dash), label ="Continuous")
plot!(m.kgrid, err2, color ="blue", line =(1.5, :solid), label ="Discrete")
ylims!(-15e-4, 5e-4)
ylabel!("Euler-equation errors")
xlabel!(L"\mathrm{Amount\ of\ capital :  k }")
#savefig("Fig3_pndp6.pdf")
display(plt6)
