function [grid] = grid_mmv(xmin, xmax, theta, num_grid)
% Function grid_mmv
%  [grid] = grid_mmv(mink, maxk, num_grid, theta)
%
% Purpose:
%  Generate grids based on the idea by
%  Maliar, Maliar and Valli (2010):
%  "Solving the Incomplete Markets Model with Aggregate Uncertainty using the Krusell-Smith Algorithm,"
%  Journal of Economic Dynamics and Control
%
%  Record of revisions:
%     Date     Programmer  Description of change
%  ==========  ==========  =====================
%  02/22/2016  T. Yamada   Original code

    % Equation (7) in Maliar et al. (2010,JEDC)
    tmp = zeros(num_grid,1);
    for i = 1:num_grid
        tmp(i) = ((i-1)/(num_grid-1)).^theta * xmax;
    end

    % adjust to [xmin,xmax]
    grid = zeros(num_grid,1);
    grid(1) = xmin;
    for i = 2:num_grid
        grid(i) = grid(i-1) + (tmp(i)-tmp(i-1))/xmax*(xmax-xmin);
    end

return;