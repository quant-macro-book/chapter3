def crra(cons,γ):
    """
    ------------------------------
    === CRRA utility function ===
    ------------------------------
    <input>
    ・cons: current consumption level
    ・γ: relative risk aversion (or inverse of the intertemporal elasticity of substitution)
    <output>
    ・util: utility
    """
    import numpy  as np
    
    if γ != 1.0:
        util = cons**(1-γ) / (1-γ)
    else:
        util = np.log(cons)
    
    return util


def mu_crra(cons, γ):
    """
    -------------------------------------------------
    === Marginal utility of CRRA utility function ===
    -------------------------------------------------
    <input>
    ・cons: current consumption level
    ・γ: relative risk aversion (or inverse of the intertemporal elasticity of substitution)
    <output>
    ・marg_util: marginal utility
    """
    marg_util = cons**(-γ)
    return marg_util


def mu_crra_inv(marg_util, γ):
    """
    --------------------------------------------
    === Inverse of marginal utility function ===
    --------------------------------------------
    <input>
    ・marg_util: marginal utility u'(c)
    ・γ: relative risk aversion (or inverse of the intertemporal elasticity of substitution)
    <output>
    ・cons: current consumption level
    """
    cons = marg_util ** (-1/γ)
    return cons


def cobb_douglas(x,y,α,z):
    """
    ------------------------------------------------------
    === Cobb-Douglass type production/utility function ===
    ------------------------------------------------------
    <input>
    ・x: capital, consumption etc.
    ・y: labor, leisure etc.
    ・α: capital share, consumption share etc.,
    ・z(optional): multiplicative factor such as TFP/preference shock
    <output>
    ・output: multiplicative factor such as TFP/preference shock.
    """
    output = z * (x**α) * (y**(1-α))
    return output


def factor_price(KoverL, tfp, α, δ):
    """
    -----------------------------------------------------------------
    === Given K/L, return interest rate and wage based on the FOC ===
    -----------------------------------------------------------------
    <input>
    ・KoverL: K/L
    ・tfp: TFP level
    ・α: capital share
    ・δ: depreciation rate
    <output>
    ・rent: interest rate
    ・wage: wage
    """
    rent = tfp * α * (KoverL**(α-1)) - δ
    wage = tfp * (1-α) * (KoverL**(α))

    return rent, wage