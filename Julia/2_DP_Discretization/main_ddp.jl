using LinearAlgebra
using Plots
using Dierckx #オイラー残差のため
using DelimitedFiles
using LaTeXStrings

include("Utils.jl")
include("Models.jl")
include("ddp.jl")

m = Models(nk = 10001, kmax = 0.5) # パラメータを含む構造体のインスタンス

# 実行
vfcn, pfcn, dif, val_tmp = ddp(m)
@time ddp(m);

#解析解
AA = (1.0-m.β)^(-1) * (log(1.0-m.α*m.β) + ((m.α*m.β)/(1.0-m.α*m.β))*log(m.α*m.β))
BB = m.α/(1.0-m.α*m.β)
v_true = AA .+ BB*log.(m.kgrid)
p_true = m.α*m.β*(m.kgrid.^m.α);

# オイラー方程式から誤差を測定
nkk = 21
kgrid2 = collect(LinRange(m.kgrid[1], m.kgrid[end], nkk))
pfcn0 = zeros(nkk)

for i in 1:nkk
    pfcn0[i] = pfcn[500*(i-1)+1]
end

cons = kgrid2.^m.α + (1.0 - m.δ)*kgrid2 - pfcn0

LHS = Utils.mu_CRRA.(cons, m.γ)
kp = pfcn0
spl = Spline1D(kgrid2, pfcn0, k = 1, bc = "extrapolate")
kpp = spl(kp)
cons = kp.^m.α + (1 - m.δ) * kp - kpp
rent = m.α*kp.^(m.α - 1) .- m.δ
RHS = m.β*(1 .+ rent) .* Utils.mu_CRRA.(cons, m.γ)
err = RHS./LHS .- 1.0;

open("err_ddp.txt", "w") do io
    writedlm(io, err)
    end

#= プロット =#
plt1 = plot(m.kgrid, vfcn,
　　        label ="Approximated solution",
            line = (1.5, :solid), color ="blue",
            legend=:bottomright)
plot!(m.kgrid, v_true,
      label = "True solution",
      line = (1.5, :dash), color ="red")

xlabel!(L"\mathrm{Amount\ of\ capital :  k }")
ylabel!(L"\mathrm{Value\  function  :  V(k)} ")
xlims!(0, m.kgrid[end])
#savefig("Fig3_dndp1.pdf")
display(plt1)

plt2 = plot(m.kgrid, pfcn,
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
#savefig("Fig3_dndp2.pdf")
display(plt2)

x = 1:m.maxiter
plt3 = plot(x, dif[1,:], label = "", line = 1.5, color="blue")
xlims!(0,250)
ylims!(0, 3.0)
ylabel!("Value function approximation error")
xlabel!("Number of iterations")
#savefig("Fig3_dndp3.pdf")
display(plt3)

x = 1:m.maxiter
plt4 = plot(x, dif[2,:], label = "", line = 1.5, color="blue")
xlims!(0,250)
ylims!(0, 4.0)
ylabel!("Policy function approximation error")
xlabel!("Number of iterations")
#savefig("Fig3_dndp4.pdf")
display(plt4)

styles = filter((s->begin
                s in Plots.supported_styles()
            end), [:solid, :dash])
styles = reshape(styles, 1, length(styles))

x = 1:m.maxiter
plt5 = plot(x, dif', line = (1.5, styles), color =["blue" "red"],
     label =["Value function" "Policy function"])
xlims!(0,250)
ylims!(0, 0.1)
ylabel!("Approximation error")
xlabel!("Number of iterations")
#savefig("Fig3_dndp5.pdf")
display(plt5)

styles = filter((s->begin
                s in Plots.supported_styles()
            end), [:dash, :dot, :dashdot, :solid])
styles = reshape(styles, 1, length(styles))

plt6 = plot(m.kgrid,val_tmp, line=(1.5, styles),
     label = ["It =1" "It =3" "It =5" "Converged"],
     color = ["blue" "red" "lightblue" "purple"])
     xlabel!(L"\mathrm{Amount\ of\ capital :  k }")
     ylabel!(L"\mathrm{Value\  function  :  V(k)} ")
     title!("Convergence of value function")
     #savefig("Fig3_dndp6.pdf")
display(plt6)
