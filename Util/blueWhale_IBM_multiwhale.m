function [X, Y, sst, krill, s, direc, steps, turns, day_vec, ars_on, good_whale_vec] = blueWhale_IBM_multiwhale(par,sst_data,krill_data,grid_pars)

% Simulates NWHALES whales with ARS search
% Once whales exit, they cannot reenter the domain

N = par.rate*par.numDays;   % EDIT kept num timesteps fixed but change rate to adjust tmp res and numdays to adjust for rate
numWhales = par.numWhales; % Rename so I don't have to type par each time

% (x,y) values - Values of cells are the X,Y-coordinates
X_domain = grid_pars.xrange(1):grid_pars.resolution:grid_pars.xrange(2); 
X_domain = repmat(X_domain(1:end-1),grid_pars.numX,1);
Y_domain = grid_pars.yrange(1):grid_pars.resolution:grid_pars.yrange(2); 
Y_domain = repmat(Y_domain(1:end-1)',1,grid_pars.numY);

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
ars_on  = zeros(numWhales,N);   % 1 if whale conducts ARS, 0 if not
gammas_E = zeros(numWhales,N);  % environmental transition prob
gammas_K = zeros(numWhales,N);  % Krill trans prob
gammas   = zeros(numWhales,N);  % environ + krill trans prob
good_whale_vec = zeros(numWhales,N);    % Which whales are in the ocean

x_new = zeros(numWhales,1); 
y_new = zeros(numWhales,1);

state2_correction_angle = pi; % Mean turning angle for S2 is pi. If whale conducts ARS in S2, need to add pi to angle selected to counteract this

% Initialize varaibles
s(:,1)       = 1;               % Start all whales in transit
direc(:,1)   = pi/2;            % State all whale pointing due north
day_vec(:,1) = par.doy_start;   % Day of year whales start on

%% Sample initial condition from bounding box
sst_on = ones(numWhales,1);     % 1 if sst has not been selected 

X(:,1) = par.start_box(1) + (par.start_box(2) - par.start_box(1)).* rand(numWhales,1); % X coord
Y(:,1) = par.start_box(3) + (par.start_box(4) - par.start_box(3)).* rand(numWhales,1); % Y coord
    
idx       = coordinateToGridCell([X(:,1),Y(:,1)],grid_pars);            % Maps (X,Y) coording in meters to matrix (column, row). Order of idx is (column, row) = (X,Y)!
idx_cell  = sub2ind([grid_pars.numX,grid_pars.numY],idx(:,2),idx(:,1)); % Maps matrix (row,column) to cell number

sst_tmp     = sst_data(:,:,par.doy_start);
krill_tmp   = krill_data(:,:,par.doy_start);

sst_loc   = sst_tmp( idx_cell );    % SST at whale locations
sst_on    = sst_on.*isnan(sst_loc); % 1 if whale at a location with no sst (i.e. on land)


while sum(sst_on)  % Resample for whales with no SST (on land)
    
    bad_idx = find(sst_on > 0);  % Whales that are on land
    
    % resample
    X(bad_idx,1) = par.start_box(1) + (par.start_box(2) - par.start_box(1)).* rand(length(bad_idx),1);
    Y(bad_idx,1) = par.start_box(3) + (par.start_box(4) - par.start_box(3)).* rand(length(bad_idx),1);
    
    idx(bad_idx,:)     = coordinateToGridCell([X(bad_idx,1),Y(bad_idx,1)],grid_pars); % Order of idx is (column, row) = (X,Y)!
    idx_cell(bad_idx)  = sub2ind([grid_pars.numX,grid_pars.numY],idx(bad_idx,2),idx(bad_idx,1));

    % Determine if new locatio is on land
    sst_loc(bad_idx)   = sst_tmp( idx_cell(bad_idx) );
    sst_on    = sst_on.*isnan(sst_loc);
    
 
end

% Record sst and krill
sst(:,1)   = sst_loc;
krill(:,1) = krill_tmp( idx_cell );

% Sample steps and turns
steps(:,1) = gamrnd(par.mu(s(:,1)).^2./par.sigma(s(:,1)).^2, par.sigma(s(:,1)).^2./par.mu(s(:,1)));
turns(:,1) = mod(vmrand(par.kappa_mu(s(:,1)),par.kappa(s(:,1))),2*pi);

%% Generate rest of sequence

run_thresh = 2*N;  % To make sure an infinite while loop isn't created
good_whales = find(sst_on == 0); % Whales in the ocean 

for k = 2:N  % Time step

    sst_on = ones(numWhales,1);  % 1 if have not yet selected sst/krill for new whale location
    
    % Generate new location for whales that were in ocean on last time step
    x_new(good_whales) = cos(direc(good_whales,k-1) + turns(good_whales,k-1)); 
    y_new(good_whales) = sin(direc(good_whales,k-1) + turns(good_whales,k-1));
    
    X(good_whales,k) = X(good_whales,k-1) + steps(good_whales,k-1).*x_new(good_whales);
    Y(good_whales,k) = Y(good_whales,k-1) + steps(good_whales,k-1).*y_new(good_whales);
    
    % Determine if whales have stepped out of ROMS domain (different than being beached on land)
    [whale_in_bounds, good_whales] = find_whales_in_bounds(X(good_whales,k), Y(good_whales,k),grid_pars,par,good_whales);  % Function to find which whales are still in bounds after most recent time step. 
    
    % Whale_in bounds is vector of {0,1}'s of length numWhales, good whales gives indices of whales in bounds
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
    sst_tmp   = sst_data(  :,:,floor(k/par.rate) + par.doy_start );
    krill_tmp = krill_data(:,:,floor(k/par.rate) + par.doy_start );
    
    sst_loc(good_whales)   = sst_tmp( idx_cell(good_whales) );
    sst_on    = whale_in_bounds.*sst_on.*isnan(sst_loc);  % 1 corresponds to whale on land - only consider whales that are still in bounds
    
    % Resample any bad locations (beached whales)
    run_times = 0;
    while (sum(sst_on) && run_times < run_thresh)  % as long as sst hasn't been selected, and run times < threshold value (to ensure finite loop)
        run_times = run_times + 1;
        
        bad_idx = find(sst_on > 0); % Beached whales 
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
    
    % Select state transition rates based on current krill and environmental conditions
    % If want different trans probs, just need to chage this function - just be sure to return the prob of being in state 2
    [gammas_K(good_whales,k), gammas_E(good_whales,k), gammas(good_whales,k)] = krill_environ_probs_multiWhale(krill_tmp( idx_cell(good_whales) ), sst_loc(good_whales), good_whales, par); % Gives probability of whale being in state 2
    s(good_whales,k) = 1 + double( rand(length(good_whales),1) < gammas(good_whales,k)); % Select state based on transition rate
    
    % Conduct ARS
    x_star = zeros(numWhales,1); y_star = zeros(numWhales,1);
    
    % Remove whales whos search radius extends out of the domain.
    ars_whales = intersect(good_whales( idx(good_whales,1) > par.search_rad + 1), good_whales( idx(good_whales,2) > par.search_rad + 1)  );
    [x_star(ars_whales), y_star(ars_whales), ars_on(ars_whales,k)] = ARS_multiWhale(idx(ars_whales,:),sst_tmp,krill_tmp, par, X_domain, Y_domain, ars_whales);  % Location of highest prob of foraging
    
    % Update for whales that are in bounds and did ARS
    x_new(ars_on(:,k)==1) = x_star(ars_on(:,k)==1) - X(ars_on(:,k)==1,k); 
    y_new(ars_on(:,k)==1) = y_star(ars_on(:,k)==1) - Y(ars_on(:,k)==1,k);
    x_new(ars_on(:,k)==1) = x_new(ars_on(:,k)==1)./(sqrt(x_new(ars_on(:,k)==1).^2 + y_new(ars_on(:,k)==1).^2)); % Normalize (needed for the acos below)
                
    state2_correction_on = (s(good_whales,k)-1).*ars_on(good_whales,k);  % If whale in S2 conducts ARS, recenter turning angle around 0
                      
    direc(good_whales,k) = sign(y_new(good_whales)).*acos(x_new(good_whales)); % Select the direction the whales moved to get from k-1 to kth location. If the whales conduct ARS, this is updated to point toward area of highest krill
    
    % Select new step lengths and turning angles    
    idx_good_whales_s    = sub2ind(size(par.mu),good_whales,s(good_whales,k) );  % Tells cell index of par.mu, sigma, kappa, kappa_mu to select for each whale. 
    % EDIT note selects new step length for all the whales in the domain
    steps(good_whales,k) = gamrnd(par.mu(idx_good_whales_s).^2./par.sigma(idx_good_whales_s).^2, par.sigma(idx_good_whales_s).^2./par.mu(idx_good_whales_s));
    turns(good_whales,k) = mod(vmrand(par.kappa_mu(idx_good_whales_s),par.kappa(idx_good_whales_s)),2*pi) + state2_correction_on.*state2_correction_angle;
    
    %EDIT add par.max_step
    steps(steps > par.max_step) = par.max_step;
    
    day_vec(good_whales,k) = floor(k/par.rate) + par.doy_start;  % Just to keep track of the day of the year (debugging purposes)
     
end





