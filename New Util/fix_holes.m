function fix_holes%(percent)
% %     backup1_file = 'erd_GOES_sst_monthly_2019_formatted.mat';
% %     load(backup1_file, 'sst');
% %     sst_backup1_temp = sst;
%     
%     backup_file = 'erd_GOES_sst_monthly_2019_formatted.mat';
%     load(backup_file, 'sst');
%     sst_backup_temp = sst;
%     
% %     mid_file = 'erd_GOES_sst_8day_2019_formatted.mat';
% %     load(mid_file, 'sst');
% %     sst_mid_temp = sst;
%     
%     data_file = 'erd_GOES_sst_1day_2019_formatted.mat';
%     load(data_file);
%         
%     % 0.5 km res?
%     [x_0, y_0, z_0] = meshgrid(1:5.5:5.5*281,1:5.5:5.5*301,1:30:356);
% %     [x_1, y_1, z_1] = meshgrid(1:5.5:5.5*281,1:5.5:5.5*301,1:8:364);
%     [x_new, y_new, z_new] = meshgrid(1:5.5:5.5*281,1:5.5:5.5*301,1:356);
%     
%     % create new matrices
%     sst_backup = interp3(x_0, y_0, z_0, sst_backup_temp, x_new, y_new, z_new);
% %     sst_mid = interp3(x_1, y_1, z_1, sst_mid_temp, x_new, y_new, z_new);
%     
%     for i=1:numel(sst)
%         if sst(i) < 0
% %             if sst_mid(i) < 0
% %                 sst_mid(i) = sst_backup(i);
% %             end
% %             sst(i) = sst_mid(i);
%             sst(i) = sst_backup(i);
%         end
%     end 
%     
%     save('test_fix_holes_input_2.mat', 'grid_pars', 'latitude', 'longitude', 'mask', 'sst', 'time')
    n=3;
    backup_file = 'input_data_geo_3_tmp_8_yr_2008.mat';
    load(backup_file, 'sst_data', 'krill_data');
    sst_backup_temp = sst_data;
    krill_backup_temp = krill_data;
    
    data_file = ['holey_input_' num2str(percent) '_percent.mat'];
    load(data_file, 'sst_data', 'krill_data', 'grid_pars');
        
    % will be a problem if grid is unclear from satellite data
    [x_0, y_0, z_0] = meshgrid(1:n:3*361,1:n:3*361,1:8:366);
    [x_new, y_new, z_new] = meshgrid(1:n:3*361,1:n:3*361,1:366);

    % create new matrices
    sst_backup   = interp3(x_0, y_0, z_0, sst_backup_temp, x_new, y_new, z_new);
    krill_backup = interp3(x_0, y_0, z_0, krill_backup_temp, x_new, y_new, z_new);
    
    for i=1:numel(sst_data)
        if sst_data(i) < 0
            sst_data(i) = sst_backup(i);
            krill_data(i) = krill_backup(i);
        end
    end 
    save(['not_' char(data_file)], 'sst_data', 'krill_data', 'grid_pars')
end