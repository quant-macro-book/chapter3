"""
k'を1つ与えた際にベルマン方程式を返す

### インプット
`params::Params`: パラメータなどをまとめた変数
`wealth::Float64`: 今期利用可能な資産
`kprime::Float64`: 次期の資本量
`vnext`: 補間した次期の価値関数

# アウトプット 
`value::Float64`: 負値にしたベルマン方程式
"""
function BellmanEq(params, wealth , kprime, vnext)
    cons = wealth - kprime
    current_util = MyEconFcn.crra(cons, params.γ)
    value = current_util + params.β*vnext(kprime)
    value = -1*value
    return value
end
