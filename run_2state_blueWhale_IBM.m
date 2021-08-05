% Blue Whale IBM
% Original code by Stephanie Dodson (sadodson@ucdavis.edu)
% Modified by Sameerah Helal (shelal@ucdavis.edu)

function run_2state_blueWhale_IBM(n, m, steps_per_day, percent)
% Run multiwhale blue whale IBM simulation for specified year    
    %% Setup
    % default: default n = 3 (spatial resolution), m = 1 (temporal),
    % steps_per_day = 4 (steps/day)
%     percent = 70;
    % check m and n to make sure it matches the size of the dataframe)
    % set par.w to [1,2,3] for default and [1,0,1} krill only
    % for test fix holes, switch grid_pars.yrange and .xrange AND
    % grid_pars.numX and .numY AND set krill to zeros(size(sst_data))
    m_1 = 1;
    years = 2008;
    % set file names and paths
    % file name head
    if percent ~= 0
        cloud_string = ['_clouded_' num2str(m_1)];
    else
        cloud_string = ['_yr_2008'];
        m_1 = m;
    end
    file_name_head = ['geo_' num2str(n) '_tmp_' num2str(m) cloud_string];
    
    % input file names
    input_file_name = ['input_data_' file_name_head '.mat'];
    
    % output file names
    if steps_per_day == 4
        rate_char = '';
    else
        rate_char = ['_rate_' num2str(steps_per_day)];
    end
    output_file_name = ['output_data_' file_name_head rate_char '.mat'];
    
    % path names
    inpath = 'C:\Users\samih\OneDrive - University of California, Davis\Percent Input Data\';
    outpath = 'C:\Users\samih\OneDrive - University of California, Davis\Percent Output Data\';

    month_val = '05'; % Starting month: 5 = May
    addpath 'C:\Users\samih\OneDrive\Documents\Blue Whale Data Resolution Project\Util'  % Add path to the functions this code calls

    % Define system parameters and store in Matlab structure 'par'
    par.rate            = steps_per_day*m_1;  % Steps per day, timesteps (1,2,3,4)
    par.numWhales       = 2000;             % Number of whales
    par.numDays         = 240/m_1;            % Number of NOT days, timesteps
    par.search_rad      = 2;                % search radius for ARS (number of grid cells)
    par.search_thres    = 0.05;             % Conduct ARS if difference between min and max on search radius are greater than this threshold
    par.w               = [1;2;3];          % Weights: Environmental, krill, normalization factor  <--- !!! Set to [1;0;1] for SST only, [0;1;1] for krill only; [1;2;3] default
    par.start_box       = [400000, 1000000, 50000, 400000]; %[700000, 700003, 100000, 100003] min(X), max(X), min(Y), max(Y)

    % safeguard
    par.max_step = 1000*100;
    
    % set hours
    hours = 24/steps_per_day;

    % mu      = [22200, 6300];    % Rate 4 mean step lengths: s1, s2
    % sigma   = [12900, 5820];    % Rate 4 std dev step lengths: s1, s2
    % mu          = [88800, 23000];  % Rate 1 mean step lengths: s1, s2
    % sigma       = [51000, 23000];  % Rate 1 std dev step lengths: s1, s2

    mu = hours*[3700, 1050];     % Rate 4 mean step lengths: s1, s2
    sigma = hours*[2150, 970];  % Rate 4 std dev step lengths: s1, s2

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


    %% Run Simulations for each year
    first_start_time = tic;
    for k = 1:1
        par.start_date  = datestr([ month_val '-01-' num2str(years(k))]);  % mm-dd-yyyy
        ref_date        = datestr(['01-01-' num2str(2008)]);
        par.doy_start   = floor((datenum(par.start_date) - datenum(ref_date))/m_1);  % Day of the year the simulation starts on

        % load the data
        load([inpath, input_file_name], 'sst_data', 'krill_data', 'grid_pars');

%         % ADD IN EDIT
%         grid_pars.resolution = 3000;
%         grid_pars.xrange  = grid_pars.resolution.*[0;size(sst_data,1)];
%         grid_pars.yrange  = grid_pars.resolution.*[0;size(sst_data,2)];
%         grid_pars.numX    = size(sst_data,1);
%         grid_pars.numY    = size(sst_data,2);
        
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

        % save output data as mat file
        save([outpath, output_file_name],'X','Y','sst','krill','s','direc','steps', 'turns','day_vec','ars_on','par','grid_pars') 
        % NOTE: if loading files into R, don't save as '-v7.3' - (remove this option)
        

        % % draw whales on map 
        %     X_domain = grid_pars.xrange(1):grid_pars.resolution:grid_pars.xrange(2); X_domain = X_domain(1:end-1);
        %     Y_domain = grid_pars.yrange(1):grid_pars.resolution:grid_pars.yrange(2); Y_domain = Y_domain(1:end-1);
        
        %     [XX, YY] = meshgrid(X_domain, Y_domain);
        
        %     writematrix([XX, YY], 'map.csv')
        %     writematrix(krill_data(:,:,245), 'krill.csv')
        
        %     figure; pcolor(XX,YY,krill_data(:,:,max(day_vec(:)))); shading interp;
        %     hold on;
        %     for j = 1:par.numWhales
        %         plot(X(j,:),Y(j,:),'k.','linewidth',2);
        %     end
        %
        %     title([num2str(n) ' km, ' num2str(m) ' days, ' num2str(steps_per_day)  ' steps/day ' '(' num2str(years(k)) ')'] ); set(gca,'fontsize',14);
        %     hold off; drawnow;

        % % save map as png image
        %     saveas(gcf,['map_geo_' num2str(n) '_tmp_' num2str(m) '_yr_2008_rate_' num2str(steps_per_day) '.png']);
        %     saveas(gcf,['not_holey_map_' num2str(percent) '_percent.png']);
        %     saveas(gcf,'test_fix_holes_plot.png');
    end

    final_stop_time = toc(first_start_time);
    disp(['Number of Whales: ' num2str(par.numWhales) ', Final Run Time: ' num2str(final_stop_time)]);
    clear all;
end
