% Put Matlab whale location data into CSV format for R test
% close all; clear all;

% getting location vectors for Sep 1, 2008
month_name = "sep";
month_day = 245; % day of the year

% % create names: time fixed, rate optional
% geo_vals = string([6;6;9;9;12;12]);%string([3; 6; 9; 12]);
% tmp_val = string([3;8;3;8;3;8]);%"8";
% rate_val_nums = [2;2;1;1;1;1];
% rate_vals = string(rate_val_nums);%string(1:2); % comment out if no rate

inpath = 'C:\Users\samih\OneDrive - University of California, Davis\Percent Output Data\';
outpath = 'C:\Users\samih\OneDrive - University of California, Davis\Percent Output Data\Sep\';
% outpath = 'C:\Users\samih\OneDrive\Documents\Blue Whale Data Resolution Project\New Output Data\September Locs\';
% filename_head = "output_data_geo_" + geo_vals + "_tmp_" + tmp_val + "_clouded_1_rate_" + rate_vals;
% input_files = inpath + filename_head + ".mat";%"_yr_2008.mat";
% output_files = outpath + filename_head + "_sep.csv";
% 
% input_files = inpath + ["output_data_geo_3_tmp_3_clouded_1.mat"; "output_data_geo_3_tmp_8_clouded_1.mat"];
% rate_val_nums = [4,4];

input_files = dir(fullfile(inpath, "*.mat"));

% loop through all files and save positions
for i = 1:numel(input_files)
    file_name = input_files(i).name;
    load(file_name);
    
    out_file_name = strrep(file_name, ".mat", ".csv");
    out_file_name = strrep(out_file_name, "output_data_", "");
    out_file_name = strrep(out_file_name, "_yr_2008", "");
    
    if contains(out_file_name,"rate_1")
        rate_val = 1;
    elseif contains(out_file_name,"rate_2")
        rate_val = 2;
    else
        rate_val = 4;
    end
    % not dividing by (4/i) gives index out of bounds
    % size(X) for rate 4: 960, rate 2: 480, rate 1: 245 instead of 240
    k = floor(month_day/(4/rate_val));
    % find whales with positive position?
    tmpX = X(:,k);
    tmpY = Y(:,k);
    idx = find(tmpX > 0);
    tmpX = tmpX(idx);
    tmpY = tmpY(idx);
    % generate csv
    writematrix([tmpX,tmpY], outpath + "sep_" + out_file_name);
    disp(out_file_name)
end
% for h = 1:numel(geo_vals)
%     for i = 1:numel(rate_vals)
%         load(input_files(h,i));
%         % not dividing by (4/i) gives index out of bounds
%         % size(X) for rate 4: 960, rate 2: 480, rate 1: 245 instead of 240
%         k = floor(month_day/(4/i));
%         % find whales with positive position?
%         tmpX = X(:,k);
%         tmpY = Y(:,k);
%         idx = find(tmpX > 0);
%         tmpX = tmpX(idx);
%         tmpY = tmpY(idx);
%         % generate csv
%         writematrix([tmpX,tmpY], output_files(h,i));
%         disp(output_files(h,i))
%     end
% end