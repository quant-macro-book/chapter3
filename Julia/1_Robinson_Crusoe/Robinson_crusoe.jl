function Robinson_crusoe(m)
    """
    最適化(optimization)と内挿法(interpolation)をつかって
    ロビンソン・クルーソー経済を解く.

    # Arguments
    -`m::Models`: パラメータを含む構造体

    # Returns
    -`vfcn::Matrix{Float64,2}`: 価値関数
    -`pfcn::Matrix{Float64,2}`: 政策関数
    -`cfcn::Matrix{Float64,2}`: 消費関数
    """

    #= 変数を定義 =#
    pfcn = zeros(m.nk, m.TT) # 政策関数
    vfcn = similar(pfcn) # 価値関数
    cfcn = similar(pfcn) # 消費関数
    wealth = m.kgrid.^m.α # 利用可能な資産

    #= 最終期(全てを消費) =#
    # pfcn[:,end] .= 0.0 # 全て消費するので貯蓄はゼロ
    cfcn[:, end] = wealth
    vfcn[:,end] = Utils.CRRA.(cfcn[:,end], m.γ) # 消費から得られる効用

    # メインループ
    @inbounds for t in m.TT-1 : -1 : 1 # 滞在期間について後ろから解いていく

        #次期の価値関数を補間
        #vnext = Spline1D(m.kgrid, vfcn[:,t+1], k = 1, bc = "extrapolate") #線形補間
        @views vnext = Spline1D(m.kgrid, vfcn[:,t+1], k = 3, bc = "extrapolate") #スプライン補間

        @inbounds for i in 1:m.nk # 状態変数についてループを行う

            BellmanEq!(kprime) = BellmanEq(m, wealth[i], kprime, vnext)
            res = optimize(BellmanEq!, 0.0, wealth[i], GoldenSection()) #最適化

            pfcn[i, t] = res.minimizer
            vfcn[i, t] = - res.minimum # 最小値を探していたので符号を反転させる
        end

        #消費関数を計算
        @views cfcn[:, t] = m.kgrid.^m.α - pfcn[:, t]
    end

    return vfcn, pfcn, cfcn
end
