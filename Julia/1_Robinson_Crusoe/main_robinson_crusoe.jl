using LinearAlgebra
using Plots
using Dierckx
using Optim
using LaTeXStrings

include("Models.jl")
include("Utils.jl")
include("BellmanEq.jl")
include("Robinson_crusoe.jl")

m = Models()
#=実行=#
vfcn, pfcn, cfcn = Robinson_crusoe(m)
@time Robinson_crusoe(m)

# 解析的解
p_true = zeros(m.nk, m.TT)

for t in 1 : m.TT
    for i in 1 : m.nk
        p_true[i, t] = m.α*m.β*((1.0-(m.α*m.β)^(m.TT-t)) / (1.0-(m.α*m.β)^(m.TT-t+1)) )*(m.kgrid[i]^m.α)
    end
end

#= プロット =#

plt1 = plot(m.kgrid, vfcn[:, end], label ="t=10",
            legend =:bottomright,show = true)
plot!(m.kgrid, vfcn[:, 9], label ="t=9")
plot!(m.kgrid, vfcn[:, 8], label ="t=8")
plot!(m.kgrid, vfcn[:, 7], label ="t=7")
plot!(m.kgrid, vfcn[:, 1], label ="t=1")
xlabel!(L"\mathrm{Amount\ of\ capital :  k_t }")
ylabel!(L"\mathrm{Value\  function  :  V_t(k_t)} ")
#savefig("Fig3_rc1.pdf")
display(plt1)


plt2 = plot(m.kgrid, pfcn[:, end], label ="t=10",
            legend =:bottomright)
plot!(m.kgrid, pfcn[:, 9], label ="t=9")
plot!(m.kgrid, pfcn[:, 8], label ="t=8")
plot!(m.kgrid, pfcn[:, 7], label ="t=7")
plot!(m.kgrid, pfcn[:, 1], label ="t=1")
xlabel!(L"\mathrm{Amount\ of\ capital\ at\ t :  k_t }")
ylabel!(L"\mathrm{Amount\ of\ capital\ at\ t+1 :  k_{t+1} }")
#savefig("Fig3_rc2.pdf")
display(plt2)


plt3 = plot(m.kgrid, cfcn[:, end], label ="t=10",
            legend =:bottomright)
plot!(m.kgrid, cfcn[:, 9], label ="t=9")
plot!(m.kgrid, cfcn[:, 8], label ="t=8")
plot!(m.kgrid, cfcn[:, 7], label ="t=7")
plot!(m.kgrid, cfcn[:, 1], label ="t=1")
xlabel!(L"\mathrm{Amount\ of\ capital\ at\ t :  k_t }")
ylabel!(L"\mathrm{Consumption :  c_t }")
#savefig("Fig3_rc3.pdf")
display(plt3)

plt4 = plot(m.kgrid, pfcn[:,9] , label="approximation at t=9",
            legend=:bottomright)
plot!(m.kgrid, p_true[:,9], label="closed form solution at t=9")
plot!(m.kgrid, pfcn[:,1], label="approximation at t=1")
plot!(m.kgrid, p_true[:,1], label="closed form solution at t=9" )
title!("True policy function and approximated solution")
xlabel!(L"\mathrm{Amount\ of\ capital\ at\ t :  k_t }")
ylabel!(L"\mathrm{Amount\ of\ capital\ at\ t+1 :  k_{t+1} }")
#savefig("Fig3_rc4.pdf")
display(plt4)

plt5 = plot(m.kgrid, vfcn[:, end-1], label ="t=10",
            legend =:bottomright,seriestype=:scatter)
xlabel!(L"\mathrm{Amount\ of\ capital :  k}")
ylabel!(L"\mathrm{Value\  function  :  V_{T-1}(k^i)} ")
xlims!(0, 1.5)
ylims!(-3.0, 0.0)
#savefig("Fig3_data.pdf")
display(plt5)

#=以下はhow_to_interpolation.mと同等の内容 =#
#= 内挿補間 =#
nkk = 1001
kkmin = 0.0
kkmax = 1.5
kkgrid = collect(LinRange(kkmin, kkmax, nkk))

v_ln = zeros(nkk)
v_cs = zeros(nkk)
ln_itp = Spline1D(m.kgrid, vfcn[:,end-1], k = 1, bc = "extrapolate") #線形補間
cs_itp = Spline1D(m.kgrid, vfcn[:,end-1], k = 3, bc = "extrapolate") #スプライン補間
for i in 1:nkk
    v_ln[i] = ln_itp(kkgrid[i])
    v_cs[i] = cs_itp(kkgrid[i])
end

plt6 = plot(kkgrid, v_ln, label="linear interpolation",
            color = "blue", legend=:bottomright);
plot!(kkgrid, v_cs, label="cubic spline interpolation",
      color = "red", line =:dash)
xlabel!(L"\mathrm{Amount\ of\ capital :  k}")
ylabel!(L"\mathrm{Value\  function  :  V_{T-1}(k)} ")
xlims!(0, 1.5)
ylims!(-3.0, 0.0)
#savefig("Fig3_interp.pdf")
display(plt6)
