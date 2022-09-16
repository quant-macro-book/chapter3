"""
状態変数と操作変数を離散化して動的計画法(discretized DP)を解く.

### Arguments
`params::Params`: パラメータを含む構造体

### Returns
`vfcn::Vector{Float64}`: 価値関数
`pfcn::Vector{Float64}`: 政策関数
`cfcn::Vector{Float64}`: 消費関数
`dif::Vector{Float64}`: 価値関数と政策関数の繰り返し計算誤差
`val_tmp::Vector{Float64}`: 価値関数の図示の際に利用
"""
function ddp(params::Params)

    count = 1 # 価値関数の図示の際に利用

    # STEP 1(b): 価値関数・政策関数の初期値を設定

    vfcn = zeros(params.nk)
    pfcn = zeros(Int64,params.nk)
    Tvfcn = zeros(params.nk)
    Tpfcn = zeros(Int64,params.nk)
    val_tmp = zeros(params.nk, 4) # プロット用の変数
    dif = zeros(2, params.maxit); # プロット用の変数

    #STEP 3: 効用関数の組み合わせ

    # 効用関数の初期値 (消費が0以下になる組み合わせにはペナルティ)
    util = -10000.0*ones(params.nk, params.nk)

    # 利用可能な最大の資源
    wealth = params.kgrid.^params.α + (1-params.δ)*params.kgrid

    # 消費が正値になる(k,k')の組み合わせについて効用を計算
    for i = 1:params.nk, j = 1:params.nk
        cons = wealth[i] - params.kgrid[j]
        if cons > 0.0
            util[j, i] = MyEconFcn.crra(cons, params.γ)
        end
    end

    # STEP 4: 価値関数を繰り返し計算
    for it = 1:params.maxit

        vkp = zeros(params.nk, params.nk)

        for i = 1:params.nk
            vkp[:, i] = util[:, i] + params.β*vfcn
        end

        for i = 1:params.nk
            Tvfcn[i] = maximum(vkp[:, i])
            Tpfcn[i] = argmax(vkp[:, i])
        end

        # 繰り返し計算誤差を確認
        dif1 = maximum(abs.((Tvfcn-vfcn)./vfcn))
        dif2 = maximum(abs.((Tpfcn-pfcn)./pfcn))

        # 価値関数・政策関数をアップデート
        vfcn = deepcopy(Tvfcn)
        pfcn = deepcopy(Tpfcn)

        # 収束途中の繰り返し計算誤差を保存
        #途中経過を図示する目的なので、通常は不要(むしろ遅くなるので消すべき)
        dif[1, it] = dif1
        dif[2, it] = dif2

        # 同じく価値関数の収束を図示する目的で保存(本来は不要)
        if it==1 || it==3 || it==5
            val_tmp[:, count] = vfcn
            count = count + 1
        end

        # println("iteration counter: $it")
        # println("error (policy): $dif1")
        # println("error (value):  $dif2")
        # println()

        if dif1 < params.tol
            val_tmp[:, end] = vfcn
            break
        end

        if it == params.maxit
            println("The model does not converge")
        end
    end

    return vfcn, pfcn, dif, val_tmp

end