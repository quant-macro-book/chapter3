function [grid] = grid_exp1(mink, maxk, num_grid)
% Function grid_exp1
%  [grid] = grid_exp1(mink, maxk, num_grid)
%
% Purpose:
%  Generate exponentially-spaced grids.
%
%  Record of revisions:
%     Date      Programmer   Description of change
%  ==========   ==========   =====================
%  02/20/2013   T. Yamada    Original code

    maxd = log(maxk+1.0);
    mesh = linspace( mink, maxd, num_grid);
    grid = exp(mesh)-1.0;

return;