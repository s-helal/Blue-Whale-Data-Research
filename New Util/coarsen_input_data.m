function coarsen_input_data(n,m)
% Each element of the array will cover from 3x3x1 to nxnxm metres^2*days

    % One year: get just one file
    %nc_file_names = dir('Data/wc15n_srf_2008.nc');  

    % From run_2state_blueWhale_IBM.m
    %nc_file = [nc_file_names.folder '/' nc_file_names.name];

    % Read data from ncdf files
    %ncread(nc_file,'temp'); sst_data_temp = permute(sst_data_temp,[2,1,4,3]);
    %ncread(nc_file,'Pzooplankton'); krill_data_temp = permute(krill_data_temp,[2,1,4,3]);
    inpath = 'C:\Users\samih\OneDrive\Documents\Blue Whale Data Resolution Project\Input Data-selected\';
    load([inpath, 'input_data_geo_3_tmp_1_yr_2008.mat'], 'sst_data', 'krill_data');    

    sst_data_temp   = sst_data;
    krill_data_temp = krill_data;
    
    % mesh the grids ORIGINAL
%     [x_0,y_0,z_0]       = meshgrid(1:3:3*361,1:3:3*361,1:366);    % 361 = size(sst_data_temp,1), etc
%     [x_new,y_new,z_new]	= meshgrid(1:n:3*361,1:n:3*361,1:m:366);

    % SPACE ONLY
    [x_0,y_0,z_0]       = meshgrid(1:3:3*361,1:3:3*361,1:366);    % 361 = size(sst_data_temp,1), etc
    [x_new,y_new,z_new]	= meshgrid(1:n:3*361,1:n:3*361,1:366);

    % create new matrices
    sst_data    = interp3(x_0,y_0,z_0,sst_data_temp,x_new,y_new,z_new);
    krill_data  = interp3(x_0,y_0,z_0,krill_data_temp,x_new,y_new,z_new);
    
%     % TIME ONLY
    time_idx = m:m:size(sst_data,3);
    sst_avg = zeros(size(sst_data,1), size(sst_data,2),length(time_idx));
    krill_avg = sst_avg;
    
    % EDIT
    for j = 1:length(time_idx)
        test_sst = sst_data(:,:,time_idx(j) - m + 1:time_idx(j));
        sst_avg(:,:,j) = mean(test_sst,3);
        test_krill = krill_data(:,:,time_idx(j) - m + 1:time_idx(j));
        krill_avg(:,:,j) = mean(test_krill,3);
    end

    sst_data = sst_avg;
    krill_data = krill_avg;
    % Setting up resolution information - corresponds the the input data
    % from the ncdf files
    grid_pars.resolution    = 1000*n; 
    grid_pars.xrange    = grid_pars.resolution.*[0;size(sst_data,1)];
    grid_pars.yrange    = grid_pars.resolution.*[0;size(sst_data,2)];
    grid_pars.numX      = size(sst_data,1);
    grid_pars.numY      = size(sst_data,2);

    outpath = 'C:\Users\samih\OneDrive\Documents\Blue Whale Data Resolution Project\New Input Data\';
    mat_file_name = 'input_data_geo_' + string(n) + '_tmp_' + string(m) + '_yr_2008.mat';

    % Save the file with variables sst_data, krill_data, grid_pars
    save(outpath+ mat_file_name, 'sst_data', 'krill_data', 'grid_pars');
    
end
