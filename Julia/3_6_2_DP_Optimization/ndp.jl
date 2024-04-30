"""
状態変数のみ離散化して操作変数は連続的に値を取る場合の動的計画法(parametric DP)の解法.
アルゴリズムの詳細は、Johnson et al. (1993)を参照

### Inputs
`params::Params`: パラメータ等を含む構造体

### Outputs
`vfcn0::Vector{Float64}`: 計算によって得られた価値関数
`pfcn1::Vector{Float64}`: 計算によって得られた政策関数
"""
function ndp(params)

    # 価値関数と政策関数の初期化
    pfcn0 = zeros(params.nk)
    vfcn0 = MyEconFcn.crra.(params.kgrid.^params.α + (1-params.δ)*params.kgrid, params.γ)

    pfcn1 = zeros(params.nk)
    vfcn1 = zeros(params.nk)

    # 利用可能な資産をあらかじめ計算しておく
    wealth = params.kgrid.^params.α + (1-params.δ)*params.kgrid

    # 繰り返し誤差を保存する変数を設定
    # 途中経過を図示する目的なので、通常は不要(むしろ遅くなるので消すべき)
    dif = zeros(2, params.maxit)

    # 価値関数を繰り返し計算
    for it = 1:params.maxit

        #次期の価値関数を補間
        #vnext = Spline1D(params.kgrid, vfcn0, k=1, bc="extrapolate") #線形補間
        vnext = Spline1D(params.kgrid, vfcn0, k=3, bc="extrapolate") #スプライン補間

        for i = 1:params.nk
            BellmanEq!(kprime) = BellmanEq(params, wealth[i], kprime, vnext)
            res = optimize(BellmanEq!, 0.0, wealth[i], GoldenSection()) # 最適化
            pfcn1[i] = res.minimizer
            vfcn1[i] = -res.minimum # 最小値を探していたので符号を反転させる
        end

        dif1 = maximum(abs.((vfcn1 - vfcn0)./vfcn0)) # 価値関数の繰り返し計算誤差
        dif2 = maximum(abs.((pfcn1 - pfcn0)./pfcn0)) # 政策関数の繰り返し計算誤差(図示のため)

        # 収束途中の繰り返し計算誤差を保存
        dif[1, it] = dif1
        dif[2, it] = dif2

        vfcn0 = deepcopy(vfcn1)
        pfcn0 = deepcopy(pfcn1)
        if dif1 < params.tol
            break
        end

        if it == params.maxit
            println("The model does not converge...")
        end
    end

    return vfcn0, pfcn0, dif
end
