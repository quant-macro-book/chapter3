module Utils

function grid_exp1(mink, maxk, num_grid)
    """
    Generate exponentially-spaced grids

    # Arguments
    - `mink::Real`: minimum value of grid points
    - `maxk::Real`: maximum value of grid points
    - `num_grid::Int`: Number of grid points

    # Return 
    - `grid::Vector`: exponentially-spaced grids
    """

    maxd = log(maxk + 1.0)
    mesh = collect(LinRange(mink, maxd, num_grid))
    grid = exp.(mesh) .- 1.0
    return grid
end

function grid_exp2(mink, maxk, num_grid)
    """
    Generate double exponentially-spaced grids

    # Arguments
    - `mink::Real`: minimum value of grid points
    - `maxk::Real`: maximum value of grid points
    - `num_grid::Int`: Number of grid points

    # Return 
    - `grid::Vector`: double exponentially-spaced grids
    """
    maxd = log(log(maxk+1)+1)
    mesh = collect(LinRange(mink, maxd, num_grid))
    grid = exp.(exp.(mesh) .- 1.0) .- 1.0

    return grid 
end

function grid_exp3(mink, maxk, num_grid)
    """
    Generate triple exponentially-spaced grids

    # Arguments
    - `mink::Real`: minimum value of grid points
    - `maxk::Real`: maximum value of grid points
    - `num_grid::Int`: Number of grid points

    # Return 
    - `grid::Vector`: triple exponentially-spaced grids
    """
    maxd = log(log(log(maxk+1)+1)+1)
    mesh = collect(LinRange(mink, maxd, num_grid))
    grid = exp.(exp.(exp.(mesh) .- 1.0) .- 1.0) .- 1.0

    return grid 
end

function grid_mmv(xmin, xmax, num_grid, theta)
    """
    Generate grid based on Maliar, Maliar and Valli(2010)
    "Solving the Incomplete Markets Model with Aggregate Uncertainty 
     using the Krusell-Smith Algorithm"
    Journal of Economic Dynamics and Control

    NB:
    In the paper, they suggest theta = 7.
    If theta = 1, grid points that are distributed uniformaly.

    # Arguments
    - `xmin::Real`: minimum value of grid points
    - `xmax::Real`: maximum value of grid points
    - `num_grid::Int`: Number of grid points
    - `theta::Real`: parameter that determines the concentration of grid points 
                     at the beginning of the interval increases.
    
    # Return
    - `grid::Vector`: constructed grid
    """
    
    # Equation (7) in Maliar et al.(2010)
    tmp = zeros(num_grid)
    for i in 1:num_grid
        tmp[i] = ((i - 1)/(num_grid - 1)).^theta * xmax
    end
    
    # adjust to [xmin, xmax]
    grid = zeros(num_grid)
    grid[1] = xmin
    
    for i in 2:num_grid
        grid[i] = grid[i-1] + (tmp[i] - tmp[i-1])/xmax *(xmax - xmin)
    end
    return grid 
end



function CRRA(cons::Real, gamma)
    """
    Compute CRRA utility function
    
    # Arguments

    - `cons::Real`: consumption value
    - `gamma::Real`: relative risk aversion
    
    # Return 
    - `util::Real`: utility value 
    """
    if gamma != 1.0
        util = cons^(1.0 - gamma) / (1.0 - gamma)
    else
        util = log(cons) 
    end
    return util
end

function mu_CRRA(cons::Real, gamma)
    """
    Compute marginal utility of CRRA-type function
    
    # Arguments 
    - "cons::VecOrMat": consumption value
    - "gamma::Real": relative risk aversion
    
    # Return
    - "mu::Real": marginal utility 
    """
    
    mu = cons^-gamma
    return mu
end



end