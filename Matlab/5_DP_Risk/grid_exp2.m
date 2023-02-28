function [grid] = grid_exp2(mink, maxk, num_grid)
% Function grid_exp2
%  [grid] = grid_exp2(mink, maxk, num_grid)
%
% Purpose:
%  Generate double exponentially-spaced grids.
%
%  Record of revisions:
%     Date      Programmer   Description of change
%  ==========   ==========   =====================
%  02/20/2013   T. Yamada    Original code

    maxd = log(log(maxk+1)+1);
    mesh = linspace( mink, maxd, num_grid);
    grid = exp(exp(mesh)-1.0)-1.0;

return;