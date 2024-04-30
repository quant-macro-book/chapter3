function [grid] = grid_exp3(mink, maxk, num_grid)
% Function grid_exp3
%  [grid] = grid_exp3(mink, maxk, num_grid)
%
% Purpose:
%  Generate triple exponentially-spaced grids.
%
%  Record of revisions:
%     Date      Programmer   Description of change
%  ==========   ==========   =====================
%  02/20/2013   T. Yamada    Original code

    maxd = log(log(log(maxk+1)+1)+1);
    mesh = linspace( mink, maxd, num_grid);
    grid = exp(exp(exp(mesh)-1.0)-1.0)-1.0;

return;