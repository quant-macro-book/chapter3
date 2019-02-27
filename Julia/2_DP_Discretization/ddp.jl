function ddp(m)
    """
    状態変数と操作変数を離散化して動的計画法(discretized DP)を解く.

    # Arguments
    `m::Models`: パラメータを含む構造体
    # Returns
    `vfcn::Vector{Float64}`: 価値関数
    `pfcn::Vector{Float64}`: 政策関数
    `cfcn::Vector{Float64}`: 消費関数
    `dif::Vector{Float64}`: 価値関数と政策関数の繰り返し計算誤差
    `val_tmp::Vector{Float64}`: 価値関数の図示の際に利用
    """

    count = 1 #価値関数の図示の際に利用

    #STEP 1(b): 価値関数・政策関数の初期値を設定

    vfcn = zeros(m.nk)
    pfcn = similar(vfcn)
    Tvfcn = similar(vfcn)
    Tpfcn = similar(vfcn)
    val_tmp = zeros(m.nk, 4)
    dif = zeros(2, 1000)

    #STEP 3: 効用関数の組み合わせ

    # 効用関数の初期値 (消費が0以下になる組み合わせにはペナルティ)
    util = -10000.0*ones(m.nk, m.nk)

    # MATLAB版と効用行列の構成が異なることに注意
    wealth = m.kgrid.^m.α + (1.0-m.δ)*m.kgrid
    # 消費が正値になる(k,k')の組み合わせについて効用を計算
    @inbounds for i in 1 : m.nk # あらゆる状態変数kについて
        @inbounds for j in 1 : m.nk # あらゆる操作変数k'について
            cons = wealth[i] - m.kgrid[j]
            if cons > 0.0
                util[i,j] = Utils.CRRA(cons, m.γ)
            end
        end
    end

    # STEP 4: 価値関数を繰り返し計算
    #　最大値を求めるコーディングがMATLAB版とは違うことに注意
    for iter in 1:m.maxiter

        @inbounds for i in 1:m.nk # あらゆる状態変数kについて
            vmin = -10000.0
            @inbounds for j in 1:m.nk # あらゆる操作変数k'について
                temp = util[i,j] + m.β*vfcn[j]
                if temp > vmin
                    Tvfcn[i] = temp
                    Tpfcn[i] = m.kgrid[j]
                    vmin = temp
                    #else
                    #break #凹性の利用
                end
            end
        end

        # 繰り返し計算誤差を確認
        dif1 = maximum(abs.((Tvfcn-vfcn)./vfcn))
        dif2 = maximum(abs.((Tpfcn-pfcn)./pfcn))

        # 価値関数・政策関数をアップデート
        vfcn = copy(Tvfcn)
        pfcn = copy(Tpfcn)

        # 収束途中の繰り返し計算誤差を保存
        #途中経過を図示する目的なので、通常は不要(むしろ遅くなるので消すべき)
        dif[1, iter] = dif1
        dif[2, iter] = dif2

        # 同じく価値関数の収束を図示する目的で保存(本来は不要)
        if iter==1 || iter==3 || iter==5
            val_tmp[:, count] = vfcn
            count = count + 1
        end

        if dif1 < m.tol
            val_tmp[:, end] = vfcn
            break
        end

        if iter == m.maxiter
            println("The model does not converge")
        end
    end
    return vfcn, pfcn, dif, val_tmp
end