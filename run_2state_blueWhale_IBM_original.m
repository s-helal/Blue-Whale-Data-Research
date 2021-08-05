% Run multiwhale blue whale IBM simulation for specified year
% Stephanie Dodson: stephanie_dodson@brown.edu
% August 31, 2018

%% Setup

close all; clear all;

outpath = 'C:\Users\samih\OneDrive - University of California, Davis\Percent Input Data\Temp\';  % Path to directory in which to save file(s) %EDIT
% mkdir(outpath);         % Makes directory (if does not already exist) % EDIT

% One year: get just one file
% nc_file_names = dir('Data/wc15n_srf_2008.nc');  % EDIT
file_names = [outpath 'input_data_geo_6_tmp_8_percent_50.mat'];
years = 2008;

month_val = '05'; % Starting month: 5 = May
addpath ../Util/  % Add path to the functions this code calls

% Define system parameters and store in Matlab structure 'par'
par.rate        =  4;           % Samples/day
par.numWhales   = 2000;         % Number of whales
par.numDays     = 240;          % Number of days
par.search_rad  = 2;            % search radius for ARS (number of grid cells)
par.search_thres = 0.05;        % Conduct ARS if difference between min and max on search radius are greater than this threshold
par.w           = [1;2;3];      % Weights: Environmental, krill, normalization factor  <--- !!! Set to [1;0;1] for SST only, [0;1;1] for krill only
par.start_box   = [400000, 1000000, 50000, 400000]; % min(X), max(X), min(Y), max(Y)

mu = [22200, 6300];     % Rate 4 mean step lengths: s1, s2
sigma = [12900, 5820];  % Rate 4 std dev step lengths: s1, s2
%mu          = [88800, 23000];  % Rate 1 mean step lengths: s1, s2
%sigma       = [51000, 23000];  % Rate 1 std dev step lengths: s1, s2

kappa       = [10, 3];          % Std dev turning angle: s1, s2
kappa_mu    = [0, pi];          % Mean turning angle: s1, s2

krill_pref  = [5, 0.3];         % Krill preferences: [steepness of transtion, krill threshold]
envir_pref  = [-1,-0.2, -16];   % Environmental (SST) preferences for the SDM 

% Uncomment/adjust noise for heterogeneous whales
par.krill_pref = repmat(krill_pref,par.numWhales,1);% + normrnd(0,0.05, par.numWhales,size(krill_pref,2)); 
par.envir_pref = repmat(envir_pref,par.numWhales,1);% + normrnd(0,0.05, par.numWhales,size(envir_pref,2));
par.kappa      = repmat(kappa,par.numWhales,1)     ;% + normrnd(0,0.05, par.numWhales,size(kappa,2));
par.kappa_mu   = repmat(kappa_mu,par.numWhales,1)  ;% + normrnd(0,0.05, par.numWhales,size(kappa_mu,2));
par.mu         = repmat(mu,par.numWhales,1)        ;% + normrnd(0,0.05, par.numWhales,size(mu,2));
par.sigma      = repmat(sigma,par.numWhales,1)     ;% + normrnd(0,0.05, par.numWhales,size(sigma,2));

% EDIT
% safeguard
par.max_step = 1000*100;

%% Run Simulations for each year

first_start_time = tic;
for k = 1:1 % EDIT
    
    disp(['Year: ' num2str(years(k))])
    
    % EDIT
%     outpath_tmp = [outpath 'year_' num2str(years(k)) '_whales_' num2str(par.numWhales) '.mat'];
    outpath_tmp = [outpath 'output_data_geo_6_tmp_8_percent_50.mat'];
    par.start_date  = datestr([ month_val '-01-' num2str(years(k))]);  % mm-dd-yyyy
    % EDIT
%     ref_date        = datestr(['01-01-' num2str(year(par.start_date))]); 
    ref_date        = datestr(['01-01-2008']); 
    par.doy_start   = datenum(par.start_date) - datenum(ref_date);  % Day of the year the simulation starts on

    % EDIT
%     nc_file = [nc_file_names(k).folder '/' nc_file_names(k).name];
    
    % Read data from ncdf files % EDIT
%     sst_data   = ncread(nc_file,'temp'); sst_data = permute(sst_data,[2,1,4,3]);
%     krill_data = ncread(nc_file,'Pzooplankton'); krill_data = permute(krill_data,[2,1,4,3]);
    load(file_names);
    % Setting up resolution information - corresponds the the input data
    % from the ncdf files
    grid_pars.resolution = 3000; 
    grid_pars.xrange  = grid_pars.resolution.*[0;size(sst_data,1)];
    grid_pars.yrange  = grid_pars.resolution.*[0;size(sst_data,2)];
    grid_pars.numX    = size(sst_data,1);
    grid_pars.numY    = size(sst_data,2);

    startTime = tic;
    [X, Y, sst, krill, s, direc, steps, turns, day_vec, ars_on, good_whale_vec] = blueWhale_IBM_multiwhale(par,sst_data,krill_data,grid_pars); % Function that runs the code
    stopTime = toc(startTime);
    
    % Quick post-processing: remove the last location of the whales that went out of bounds because this location is out of bounds 
    for j = 1:par.numWhales
       
        if nnz(X(j,:)) < size(X,2)  % If there is at least one zero
            X(j,nnz(X(j,:))) = 0;
            Y(j,nnz(Y(j,:))) = 0;
        end
        
    end
   
    disp(['Run Time: ' num2str(stopTime)]);
    
    save(outpath_tmp ,'X','Y','sst','krill','s','direc','steps', 'turns','day_vec','ars_on','par','grid_pars','-v7.3'); % Save a Matlab file. NOTE: if loading files into R, don't save as '-v7.3' - (remove this option)
% EDIT
%     % Plot all the locations of the whales 
%     X_domain = grid_pars.xrange(1):grid_pars.resolution:grid_pars.xrange(2); X_domain = X_domain(1:end-1);
%     Y_domain = grid_pars.yrange(1):grid_pars.resolution:grid_pars.yrange(2); Y_domain = Y_domain(1:end-1);
% 
%     [XX, YY] = meshgrid(X_domain, Y_domain);
%     figure; pcolor(XX,YY,krill_data(:,:,max(day_vec(:)))); shading interp;
%     hold on;
%     for j = 1:par.numWhales
%         plot(X(j,:),Y(j,:),'k.','linewidth',2);
%     end
%     
%     title(['Year: ' num2str(years(k))] ); set(gca,'fontsize',18);
%     hold off; drawnow;
%    
    %saveas(gcf,['all_100whales_rho_0p3_sr8_rate1_year_' num2str(years(k))'.png']); 
    
end


final_stop_time = toc(first_start_time);

disp(['Number of Whales: ' num2str(par.numWhales) ', Final Run Time: ' num2str(final_stop_time)]);













