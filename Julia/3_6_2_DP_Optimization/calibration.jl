function calibration()
    β = 0.96
    γ = 1.0
    α = 0.4
    δ = 0.08

    nk = 101
    kmax = 10.0
    kmin = 0.05

    # 自作のコードで等分のグリッドを計算
    kgrid = GenerateGrid.grid_uni(kmin, kmax, nk)
    # これまで通り⬇でもOK
    #kgrid = collect(LinRange(kmin, kmax, nk))

    maxit = 1000
    tol = 1e-5

    return Params(β, γ, α, δ, nk, kmax, kmin, kgrid, maxit, tol)
end
