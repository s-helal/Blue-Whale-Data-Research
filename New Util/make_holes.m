function make_holes
%     data_file = ['input_data_geo_' num2str(n) '_tmp_1_yr_2008.mat'];
    data_file = 'input_data_geo_3_tmp_1_yr_2008.mat';
    load(data_file, 'sst_data', 'krill_data', 'grid_pars');
    
    elements = numel(sst_data);
    holes = floor(0.3*elements);
    locs = randi(elements, holes, 1);
    
    sst_data(locs) = -1;
    krill_data(locs) = -1;   
    
    save(['holey_input_' num2str(percent) '_percent.mat'], 'sst_data', 'krill_data', 'grid_pars')
end