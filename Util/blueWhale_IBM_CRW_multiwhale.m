function [X, Y, sst, krill, s, direc, steps, turns, day_vec, good_whale_vec] = blueWhale_IBM_CRW_multiwhale(X0,Y0,par,sst_data,krill_data,grid_pars)
% Simulates par.numWhales whale with CRW not based on krill or SST data
% Started from (X0,Y0) locations
% Once whales exit, they cannot reenter the domain
% 2 state model: transit and foraging

N = par.rate*par.numDays;
numWhales = par.numWhales;

% Preallocate matrices: Rows are whales, columns are time steps
X       = zeros(numWhales,N);   % X coordinate (in meters from south-west corner of domain)
Y       = zeros(numWhales,N);   % Y coordinate (in meters from south-west corner of domain)
sst     = zeros(numWhales,N);   % SST
krill   = zeros(numWhales,N);   % Krill
day_vec = zeros(numWhales,N);   % Day of year
s       = zeros(numWhales,N);   % State (1 or 2)
direc   = zeros(numWhales,N);   % Direction whale moved to get to current location
steps   = zeros(numWhales,N);   % Step lengths
turns   = zeros(numWhales,N);   % Turning angles
good_whale_vec = zeros(numWhales,N);    % Which whales are in the ocean

x_new = zeros(numWhales,1); 
y_new = zeros(numWhales,1);

% Initialize varaibles
s(:,1)     = 1;           % Start in transit
direc(:,1) = pi/2;    % Due north
day_vec(:,1)= par.doy_start;

% Initial locations: Assumed to be in the ocean
X(:,1) = X0;
Y(:,1) = Y0;
idx       = coordinateToGridCell([X(:,1),Y(:,1)],grid_pars);            % Maps (X,Y) coording in meters to matrix (column, row). Order of idx is (column, row) = (X,Y)!
idx_cell  = sub2ind([grid_pars.numX,grid_pars.numY],idx(:,2),idx(:,1)); % Maps matrix (row,column) to cell number

% EDIT ERROR
sst_tmp     = sst_data(:,:,min(par.doy_start,end));
krill_tmp   = krill_data(:,:,min(par.doy_start,end));

% Record sst and krill
sst(:,1)   = sst_tmp( idx_cell );
krill(:,1) = krill_tmp( idx_cell );

% Sample steps and turns
steps(:,1) = gamrnd(par.mu(s(:,1)).^2./par.sigma(s(:,1)).^2, par.sigma(s(:,1)).^2./par.mu(s(:,1)));
turns(:,1) = mod(vmrand(par.kappa_mu(s(:,1)),par.kappa(s(:,1))),2*pi);
        

%% Generate rest of sequence
run_thresh = 2*N;  % To make sure not an infinite while loop
good_whales = 1:numWhales; % Whales in the ocean 


for k = 2:N
    
    sst_on = ones(numWhales,1);
    
    % Generate new location for whales that were in ocean on last time step
    x_new(good_whales) = cos(direc(good_whales,k-1) + turns(good_whales,k-1)); 
    y_new(good_whales) = sin(direc(good_whales,k-1) + turns(good_whales,k-1));
    
    X(good_whales,k) = X(good_whales,k-1) + steps(good_whales,k-1).*x_new(good_whales);
    Y(good_whales,k) = Y(good_whales,k-1) + steps(good_whales,k-1).*y_new(good_whales);
    
    % Determine if whales have stepped out of ROMS domain (different than being beached on land)
    [whale_in_bounds, good_whales] = find_whales_in_bounds(X(good_whales,k), Y(good_whales,k),grid_pars,par,good_whales);  % Function to find which whales are still in bounds after most recent time step. 
    %Whale_in bounds is vector of {0,1}'s of length numWhales, good whales gives indices of whales in bounds
    good_whale_vec(good_whales,k) = 1; % Mostly for debugging purposes
    
     % Check to make sure there are still some whales in the domain
    if isempty(good_whales)
       disp(['All Whales Out of Bounds on Time Step ' num2str(k)]); 
       break;
        
    end
    
    % Extract SST/Krill at new location
    idx = zeros(numWhales,2); idx_cell = zeros(numWhales,1);
    sst_loc = zeros(numWhales,1);
    
    idx(good_whales,:)      = coordinateToGridCell([X(good_whales,k),Y(good_whales,k)],grid_pars); 
    idx_cell(good_whales)   = sub2ind([grid_pars.numX,grid_pars.numY],idx(good_whales,2),idx(good_whales,1));
    sst_tmp   = sst_data(  :,:,min(floor(k/par.rate) + par.doy_start ,size(sst_data,3)));
    krill_tmp = krill_data(:,:, min(floor(k/par.rate) + par.doy_start ,size(sst_data,3)));
    % EDIT ERROR above 2 used to have this: floor(k/par.rate) + par.doy_start but
    % gave error for 3rd element being too large
    sst_loc(good_whales) = sst_tmp(idx_cell(good_whales));
    sst_on = whale_in_bounds.*sst_on.*isnan(sst_loc);
    
    % Resample any bad locations (beached whales)
   run_times = 0;
   
   while (sum(sst_on) && run_times < run_thresh)
       run_times = run_times + 1;
       
       bad_idx = find(sst_on >0);  % Beached whales
       idx_bad_whales_s = sub2ind(size(par.mu),bad_idx,s(bad_idx,k-1) ); % Find the cell index
    
       % Select new steps and turns
        steps(bad_idx,k-1) = gamrnd(par.mu(idx_bad_whales_s).^2./par.sigma(idx_bad_whales_s).^2, par.sigma(idx_bad_whales_s).^2./par.mu(idx_bad_whales_s));
        %turns(bad_idx,k-1) = mod(vmrand(par.kappa_mu(idx_bad_whales_s),par.kappa(idx_bad_whales_s)),2*pi);
    
       % EDIT steps not greater than max
        steps(steps > par.max_step) = par.max_step;
        
        turns(bad_idx,k-1) = sign(-1 + 2*rand(length(bad_idx),1) )*pi/2 + turns(bad_idx,k-1); % Add or subtract 90 degrees with prob 1/2 each 
        
        % Calculate new positions
        x_new(bad_idx) = cos(direc(bad_idx,k-1) + turns(bad_idx,k-1));
        y_new(bad_idx) = sin(direc(bad_idx,k-1) + turns(bad_idx,k-1));
        X(bad_idx,k) = X(bad_idx,k-1) + steps(bad_idx,k-1).*x_new(bad_idx);
        Y(bad_idx,k) = Y(bad_idx,k-1) + steps(bad_idx,k-1).*y_new(bad_idx);
        
        idx(bad_idx,:) = coordinateToGridCell([X(bad_idx,k),Y(bad_idx,k)],grid_pars); 
        idx_cell(bad_idx)  = sub2ind([grid_pars.numX,grid_pars.numY],idx(bad_idx,2),idx(bad_idx,1));

        sst_loc(bad_idx)   = sst_tmp( idx_cell(bad_idx) );
        sst_on    = whale_in_bounds.*sst_on.*isnan(sst_loc); 
        
   end
   
    % record sst and krill
    sst(good_whales,k)   = sst_loc(good_whales);
    krill(good_whales,k) = krill_tmp( idx_cell(good_whales) );
       
   % Select state based on constant transition rates         
     s(good_whales,k) = 1 + double( rand(length(good_whales),1) < par.tau); % tau = probability of foraging
     direc(good_whales,k) = sign(y_new(good_whales)).*acos(x_new(good_whales));
             
    
   % Select new step lengths and turning angles
   idx_good_whales_s    = sub2ind(size(par.mu),good_whales,s(good_whales,k) );  % Tells cell index of par.mu, sigma, kappa, kappa_mu to select for each whale. 
   % EDIT note selects new step length for all the whales in the domain
   steps(good_whales,k) = gamrnd(par.mu(idx_good_whales_s).^2./par.sigma(idx_good_whales_s).^2, par.sigma(idx_good_whales_s).^2./par.mu(idx_good_whales_s));
   turns(good_whales,k) = mod(vmrand(par.kappa_mu(idx_good_whales_s),par.kappa(idx_good_whales_s)),2*pi);
   
   % EDIT add par.max_step
   steps(steps > par.max_step) = par.max_step;
    
   day_vec(good_whales,k) = floor(k/par.rate) + par.doy_start;  % Just to keep track of the day of the year (debugging purposes)
    
    
    
   
     
    
end






