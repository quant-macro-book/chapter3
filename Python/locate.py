def locate(xx,x):
    """
    -------------------------------------------
    === Find the position of x in vector xx ===
        ※Translation from Numerical Recipes
    -------------------------------------------
    <input>
    ・xx: vector
    ・x: points we want to know
    <output>
    ・loc: location where variable x is in vector xx.
    """
    import numpy as np

    l = len(xx)
    jl = -1
    ju = l 

    for i in range(l):
        
        if ju - jl <= 1:
            break
        jm = int(np.floor((ju+jl)/2))
        if x >= xx[jm]:
            jl = jm
        else:
            ju = jm
        
    if x == xx[1]:
        loc = 0
    elif x == xx[l-1]:
        loc = l - 2
    else:
        loc = jl

    return loc