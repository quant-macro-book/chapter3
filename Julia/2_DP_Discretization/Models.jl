struct Models{TI<:Integer, TF<:AbstractFloat, TV<:Vector}
    β::TF # 割引因子
    γ::TF # 相対的危険回避度
    α::TF # 資本分配率
    δ::TF # 固定資本減耗
    TT::TI # 無人島に滞在する期間
    nk::TI # 資本グリッドの個数
    kgrid::TV # 資本グリッド
    maxiter::TI # 繰り返し計算の最大値
    tol::TF # 許容誤差
end

function Models(;
                   β = 0.96,
                   γ = 1.0,
                   α = 0.4,
                   δ = 1.0, # 0.08
                   TT = 10,
                   nk = 11, # ロビンソンクルーソー
                        # 10001 操作変数が離散のケース
                        # 21 操作変数が連続なケース
                   kmax = 1.0,# ロビンソンクルーソー
                          # 0.5 無限期間でδ=1.0のケース
                          # 10.0 無限期間でδ=0.08のケース
                   kmin = 0.05,
                   method = "uniform" # "exp1", "exp2", "exp3", "mmv"
                    )
    
    # 資本グリッド生成
    if method == "uniform"
        kgrid = collect(LinRange(kmin, kmax, nk))
    elseif method == "exp1"
        kgrid = Utils.grid_exp1(kmin, kmax, nk)
    elseif method == "exp2"
        kgrid = Utils.grid_exp2(kmin, kmax, nk)
    elseif method == "exp3"
        kgrid = Utils.grid_exp3(kmin, kmax, nk)
    elseif method == "mmv"
        kgrid = Utils.grid_mmv(kmin, kmax, nk, 7)
    end
    
    maxiter = 1000
    tol = 1e-5
    Models(β, γ, α, δ, TT, nk, kgrid, maxiter, tol)
end