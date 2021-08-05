function create_clouds(n,m_0)%,m_1,percent_clouds,method)
    m_1 = 1;
    percent_clouds = 0.7;
    method = "";
    % n = spatial resolution, m_0 = coarse data temporal resolution; m_1 = fine data temporal resolution
    % percent_clouds (0.3,0.5,0.7,etc.) = percent of the fine data to replace with coarse data
    % method = "rep" for repetition instead of interpolation, otherwise ""
    % example create_clouds(6,8,1,0.7,"") will create 6 km with 70% 8 day
    % clouds and 30% 1 day
    
    inpath = 'C:\Users\samih\OneDrive - University of California, Davis\Input Data\';
    main_file = ['input_data_geo_' num2str(n) '_tmp_' num2str(m_1) '_yr_2008.mat'];
    load([inpath main_file], 'sst_data', 'krill_data', 'grid_pars');
    
    sst_main = sst_data;
    krill_main = krill_data;
    grid_pars_main = grid_pars;

    backup_file = ['input_data_geo_' num2str(n) '_tmp_' num2str(m_0) '_yr_2008.mat'];
    load([inpath backup_file], 'sst_data', 'krill_data');
    sst_backup_temp = sst_data;
    krill_backup_temp = krill_data;
    
    if method == "rep"
        % with repetition
        sst_backup = zeros(size(sst_main));
        krill_backup = sst_backup;

        for i = 1:size(sst_backup_temp,3)
            sst_backup(:,:,i:i+m_0-1) = repmat(sst_backup_temp(:,:,i),1,1,m_0);
            krill_backup(:,:,i:i+m_0-1) = repmat(krill_backup_temp(:,:,i),1,1,m_0);
        end
    else
        % with interpolation
        % will be a problem if grid is unclear from satellite data
        [x_0, y_0, z_0]         = meshgrid(1:n:3*361,1:n:3*361,1:m_0:size(sst_backup_temp,3)*m_0); %size(..)*m_.. was 366 for both but doesn't work with moving mean
        [x_new, y_new, z_new]   = meshgrid(1:n:3*361,1:n:3*361,1:m_1:size(sst_main,3)*m_1);

        % create new matrices
        sst_backup   = interp3(x_0, y_0, z_0, sst_backup_temp, x_new, y_new, z_new);
        krill_backup = interp3(x_0, y_0, z_0, krill_backup_temp, x_new, y_new, z_new);     
    end

    % find cloud locations
    elements = numel(sst_main);
    holes = floor(percent_clouds*elements);
    locs = randi(elements, holes, 1);
    
    sst_main(locs) = -1;
    
    % create clouds
    for i=1:numel(sst_main)
        if sst_main(i) == -1
            sst_main(i) = sst_backup(i);
            krill_main(i) = krill_backup(i);
        end
    end 
    
    % store data
    sst_data = sst_main;
    krill_data = krill_main;
    grid_pars = grid_pars_main;
    
    outpath = 'C:\Users\samih\OneDrive - University of California, Davis\Percent Input Data\';
    % add rep to name if using repetition method instead of interpolating
    % (interpolation is the default)
    if method == "rep"
        outfile = [outpath 'input_data_geo_' num2str(n) '_tmp_' num2str(m_0) '_percent_' num2str(percent_clouds*100) '_rep'];
    else
%         outfile = [outpath 'input_data_geo_' num2str(n) '_tmp_' num2str(m_0) '_percent_' num2str(percent_clouds*100)];
        outfile = [outpath 'input_data_geo_' num2str(n) '_tmp_' num2str(m_0) '_clouded_' num2str(m_1) '.mat'];
    end
    
    % default is fine data has 1 day temporal, but if not:
    if m_1 ~= 1
        % outfile = [outpath 'input_data_geo_' num2str(n) '_tmp_' num2str(m_0) '_clouded_' num2str(m_1) '_percent_' num2str(percent_clouds*100)];
        outfile = [outpath 'input_data_geo_' num2str(n) '_tmp_' num2str(m_0) '_clouded_' num2str(m_1) '.mat'];
    end
    
    save(outfile, 'sst_data', 'krill_data', 'grid_pars')
end