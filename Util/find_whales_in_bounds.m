function [whale_in_bounds, good_whales] = find_whales_in_bounds(X,Y,grid_pars,par,prev_good_whales)
% (X,Y) coordinates of whales
% Return vector of 0,1's, with 1 = whale in bounds, 0 = whale out of bounds
% Prev_good_whales: gives indicies of previously good whales

whale_in_bounds = zeros(par.numWhales,1);
whale_in_bounds(prev_good_whales) = 1;

x_left  = find(X - grid_pars.xrange(1) < 0); % < 0 if whale west of domain
x_right = find(X - grid_pars.xrange(2) > 0); % > 0 if whale east of domain

y_below = find(Y - grid_pars.yrange(1) < 0); % < 0 if whale south of domain
y_above = find(Y - grid_pars.yrange(2) > 0); % > 0 if whale north of domain

oob_idx = prev_good_whales([x_left; x_right; y_below; y_above]);  % Whales that were previously good, but have now stepped out of bounds
whale_in_bounds(oob_idx) = 0;
good_whales = find(whale_in_bounds > 0); 



