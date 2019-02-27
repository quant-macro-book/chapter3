function BellmanEq(m, wealth , kprime, vnext)
    """
    k'を1つ与えた際にベルマン方程式を返す
    robinson_crusoe.jl, main_ndp.jlから呼び出して使う
    
    # Arguments
    `m::Models`: パラメータなどを含むオブジェクト
    `wealth::Real`: 今期利用可能な資産
    `kprime::Real`: 次期の資本量
    `vnext`: 補間した次期の価値関数
    
    # Return 
    `value::Real`:　負値にしたベルマン方程式
    """
    value = - Utils.CRRA((wealth - kprime), m.γ) - m.β * vnext(kprime)
    
    return value 
end