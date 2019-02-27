function ndp(m)
    """
    状態変数のみ離散化して操作変数は連続的に値を取る場合の動的計画法(parametric DP)の解法.
    アルゴリズムの詳細は、Johnson et al. (1993)を参照

    # Arguments
    `m::Models`: パラメータ等を含む構造体

    # Returns
    `vfcn0::Vector{Float64}`: 計算によって得られた価値関数
    `pfcn1::Vector{Float64}`: 計算によって得られた政策関数
    """
    #　価値関数と政策関数の初期化

    pfcn0 = zeros(m.nk) #単純な計算には不要 プロットのためだけに用意
    vfcn0 = Utils.CRRA.(m.kgrid.^m.α + (1- m.δ)* m.kgrid, m.γ)

    pfcn1 = similar(pfcn0) # pfcn0を入れない場合は zeros(m.nk)で初期化
    vfcn1 = similar(pfcn1)

    # 利用可能な資産をあらかじめ計算しておく
    wealth =  m.kgrid.^m.α + (1 - m.δ) * m.kgrid

    # 繰り返し誤差を保存する変数を設定
    # 途中経過を図示する目的なので、通常は不要(むしろ遅くなるので消すべき)
    dif = zeros(2, m.maxiter)

    #= 価値関数を繰り返し計算 =#
    for iter in 1 : m.maxiter

        #次期の価値関数を補間
        #vnext = Spline1D(m.kgrid, vfcn0, k = 1, bc = "extrapolate") #線形補間
        vnext = Spline1D(m.kgrid, vfcn0, k = 3, bc = "extrapolate") #スプライン補間

        #temp = m.kgrid[1] #単調性を利用する際に使う
        for i in 1 : m.nk
            BellmanEq!(kprime) = BellmanEq(m, wealth[i], kprime, vnext)
            res = optimize(BellmanEq!, 0.0, wealth[i], GoldenSection()) #最適化
            #res = optimize(BellmanEq!, temp, wealth[i], GoldenSection()) #単調性の利用
            pfcn1[i] = res.minimizer
            vfcn1[i] = - res.minimum # 最小値を探していたので符号を反転させる
        end

        dif1 = maximum(abs.((vfcn1 - vfcn0)./vfcn0))　# 価値関数の繰り返し計算誤差
        dif2 = maximum(abs.((pfcn1 - pfcn0)./pfcn0)) # 政策関数の繰り返し計算誤差(図示のため)

        # 収束途中の繰り返し計算誤差を保存
        dif[1,iter] = dif1
        dif[2,iter] = dif2

        vfcn0 = copy(vfcn1)
        pfcn0 = copy(pfcn1)
        if dif1 < m.tol
            break
        end

        if iter == m.maxiter
            println("The model does not converge...")
        end
    end

    return vfcn0, pfcn0, dif
end
