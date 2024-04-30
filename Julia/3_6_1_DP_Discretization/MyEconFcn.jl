module MyEconFcn

#=
Collect functions used in Economics.
=#

export crra, mu_crra, mu_crra_inv, cobb_douglas, factor_price, accidental_bequest

"""
crra(cons::Float64, γ::Float64)

Purpose:
CRRA utility function.

Input:
cons → current consumption level,
γ → relative risk aversion (or inverst of the intertemporal elasticity of substitution).

Output:
util → utility.
"""
function crra(cons::Float64, γ::Float64)
    if γ != 1.0
        util = cons^(1-γ) / (1-γ)
    else
        util = log(cons)
    end
    return util
end


"""
mu_crra(cons::Float64, γ::Float64)

Purpose:
Marginal utility of CRRA utility function.

Input:
cons → current consumption level,
γ → relative risk aversion (or inverst of the intertemporal elasticity of substitution).

Output:
marg_util → marginal utility.
"""
function mu_crra(cons::Float64, γ::Float64)
    marg_util = cons^(-γ)
    return marg_util
end


"""
mu_crra_inv(marg_value::Float64, γ::Float64)

Purpose:
Inverse of marginal utility function.

Input:
marg_util → marginal utility u'(c),
γ → relative risk aversion (or inverst of the intertemporal elasticity of substitution).

Output:
cons → current consumption level.
"""
function mu_crra_inv(marg_util::Float64, γ::Float64)
    cons = marg_util^(-1/γ)
    return cons
end


"""
cobb_douglas(x::Float64, y::Float64, α::Float64, z::Float64=1)

Purpose:
Cobb-Douglass type production/utility function.

Input:
x → capital, consumption etc.,
y → labor, leisure etc.,
α → capital share, consumption share etc.,
z (optional) → multiplicative factor such as TFP/preference shock.

Output:
output → current output/utility level.
"""
function cobb_douglas(x::Float64, y::Float64, α::Float64, z::Float64=1.0)
    output = z*(x^α)*(y^(1-α))
    return output
end


"""
factor_price(KoverL::Float64, tfp::Float64, α::Float64, δ::Float64)

Purpose:
Given K/L, return interest rate and wage based on the FOC conditions.

Input:
KoverL → K/L,
tfp → TFP level,
α → capital share,
δ → depreciation rate,

Output:
rent → interest rate
wage → wage
"""
function factor_price(KoverL::Float64, tfp::Float64, α::Float64, δ::Float64)
    rent = tfp*   α *(KoverL^(α-1)) - δ
    wage = tfp*(1-α)*(KoverL^ α)
    return rent, wage
end

end # module
