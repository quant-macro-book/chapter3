"""
locate(xx::Array{Float64}, x::Float64)

Purpose:
Find the position of x in vector xx.

Translation from Numerical Recipes.

Input:
xx → vector,
x → points we want to know.

Output:
loc → location where variable x is in vector xx.
"""
function locate(xx::Array{Float64}, x::Float64)

    l = length(xx)
    ascnd = xx[l] >= xx[1]
    jl = 0
    ju = l + 1

    for i = 1:l
        if ju - jl <= 1
            break
        end
        jm = Int(floor((ju+jl)/2))
        if x >= xx[jm]
            jl = jm
        else
            ju = jm
        end
    end

    if x == xx[1]
        loc = 1
    elseif x == xx[l]
        loc = l - 1
    else
        loc = jl
    end

    return loc

end
