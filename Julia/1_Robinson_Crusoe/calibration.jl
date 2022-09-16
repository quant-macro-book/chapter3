function calibration()
    β = 0.96
    γ = 1.0
    α = 0.4
    δ = 1.0
    T = 10

    nk = 11
    kmax = 1.0
    kmin = 0.05

    # 自作のコードで等分のグリッドを計算
    kgrid = GenerateGrid.grid_uni(kmin, kmax, nk)
    # これまで通り⬇でもOK
    #kgrid = collect(LinRange(kmin, kmax, nk))

    return Params(β, γ, α, δ, T, nk, kmax, kmin, kgrid)
end
