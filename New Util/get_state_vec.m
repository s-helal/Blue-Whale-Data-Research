function get_state_vec%(n,m,m_0)
%% get percent foraging for clouded-- run alone
    % with clouds n=3,m=1; n=3,m=3; n=12,m=1; n=12,m=3
    % without clouds all above + n=3,m=8; n=12; m=8
    % n=3;m=1;
    % NOTE: m_0 is the coarse temporal resolution of the clouding
    
%     percent_clouds = 30;
%     inpath = "C:\Users\samih\OneDrive\Documents\Blue Whale Data Resolution Project\New Clouded Data\";
%     name_id_cloud = [num2str(n) "_km_" num2str(m) "_" num2str(m_0) "_day_" num2str(percent_clouds) "percent"];
%     cloud_file = [inpath "output_" name_id_cloud ".mat"];
%     load(cloud_file)
%     
%     for k = 1:size(Y,1)*size(Y,2)
%         if isnan(Y(k))
%         Y(k) = 0;
%         end
%     end
%     state = zeros(size(s,2),1);
%     for k=1:size(Y,2)
%         state(k) = size(nonzeros(s(:,k) == 2),1)/size(nonzeros(s(:,k)),1);
%     end
%     writematrix(state, [inpath "state_" name_id_cloud ".csv"]);
    
    inpath = "C:\Users\samih\OneDrive - University of California, Davis\Percent Output Data\";
    outpath = "C:\Users\samih\OneDrive - University of California, Davis\Percent Output Data\State\";
%     name_id_cloud = ["geo_" num2str(n) "_tmp_" num2str(m) "_clouded_" num2str(m_0)];
%     cloud_file = [inpath "output_data_" name_id_cloud ".mat"];
%     load(cloud_file)
    
    files = dir(fullfile(inpath, "*.mat"));
%     files = ["output_data_geo_6_tmp_8_percent_70.mat";
%         "output_data_geo_6_tmp_8_percent_50.mat";
%         "output_data_geo_6_tmp_8_percent_30.mat"];%; "output_data_geo_12_tmp_1_rate_1_yr_2008.mat"; "output_data_geo_6_tmp_1_rate_2_yr_2008.mat"];
    for i=1:size(files,1)
        file_name = files(i).name;
        name_id = strrep(file_name, ".mat", ".csv");
        name_id = strrep(name_id, "output_data_", "");
        name_id = strrep(name_id, "_yr_2008", "");
        load(inpath + file_name);
        for k = 1:size(Y,1)*size(Y,2)
            if isnan(Y(k))
            Y(k) = 0;
            end
        end
        state = zeros(size(s,2),1);
        for k=1:size(Y,2)
            state(k) = size(nonzeros(s(:,k) == 2),1)/size(nonzeros(s(:,k)),1);
        end
        writematrix(state, [outpath + "state_" + name_id]);
%         plot(state);
    end
    
%     inpath2 = "C:\Users\samih\OneDrive\Documents\Blue Whale Data Resolution Project\";
%     name_id_reg = ["geo_" num2str(n) "_tmp_" num2str(m) "_yr_2008"];
%     reg_file = [inpath2 "output_data_" name_id_reg ".mat"];
%     load(reg_file)
%     
%     for k = 1:size(Y,1)*size(Y,2)
%         if isnan(Y(k))
%         Y(k) = 0;
%         end
%     end
%     state = zeros(size(s,2),1);
%     for k=1:size(Y,2)
%         state(k) = size(nonzeros(s(:,k) == 2),1)/size(nonzeros(s(:,k)),1);
%     end
%     writematrix(state, [inpath "state_" name_id_reg ".csv"]);
end