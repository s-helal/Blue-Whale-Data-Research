function cellNum = coordinateToGridCell(X,grid_pars)

% Take in vector X of spatial coordinates (in meters) and find corresponding 
% cell number in ROMS data set
% Somewhat stupidly hard coded for ROMS data
% Gives (x,y) = (col,row)

cellNum = zeros(size(X));

xaxis = grid_pars.xrange(1) : grid_pars.resolution : grid_pars.xrange(2); xaxis = xaxis(1:end-1);
yaxis = grid_pars.yrange(1) : grid_pars.resolution : grid_pars.yrange(2); yaxis = yaxis(1:end-1);

tmp1 = X(:,1) - xaxis; tmp1(tmp1 > 0) = 1; tmp1(tmp1 <=0 ) = 0;
tmp2 = X(:,2) - yaxis; tmp2(tmp2 > 0) = 1; tmp2(tmp2 <=0 ) = 0;

cellNum(:,1) = sum(tmp1,2);  % Sum the number of non-zeros -- gives index number
cellNum(:,2) = sum(tmp2,2);
    
% If position is out of bounds, will give a 0, so give a grid cell number that will give a NaN value. (on land) 
cellNum(cellNum(:,1)==0,1) = length(xaxis); 
cellNum(cellNum(:,2)==0,2) = length(yaxis);  

















