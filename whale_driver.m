inpath = 'C:\Users\samih\OneDrive - University of California, Davis\Percent Input Data\Temp\';
files = dir(fullfile(inpath, "*.mat"));
for i = 1:length(files)
    cur_file = files(i).name;
    disp(cur_file);
    make_map(cur_file)
    close all;
%     figure; pcolor(XX,YY,krill_data(:,:,245)); shading interp;
end
% generate output data
% run_2state_blueWhale_IBM(3,1,4,0);
% run_2state_blueWhale_IBM(3,1,2,0);
% run_2state_blueWhale_IBM(3,1,1,0);
% run_2state_blueWhale_IBM(6,1,4,0);
% run_2state_blueWhale_IBM(6,1,2,0);
% run_2state_blueWhale_IBM(6,1,1,0);
% run_2state_blueWhale_IBM(9,1,4,0);
% run_2state_blueWhale_IBM(9,1,2,0);
% run_2state_blueWhale_IBM(9,1,1,0);
% run_2state_blueWhale_IBM(12,1,4,0);
% run_2state_blueWhale_IBM(12,1,2,0);
% run_2state_blueWhale_IBM(12,1,1,0);

% run_2state_blueWhale_IBM(3,3,4,0);
% run_2state_blueWhale_IBM(3,3,2,0);
% run_2state_blueWhale_IBM(3,3,1,0);
% run_2state_blueWhale_IBM(3,3,4,70);
% run_2state_blueWhale_IBM(3,3,2,70);
% run_2state_blueWhale_IBM(3,3,1,70);
% run_2state_blueWhale_IBM(6,3,4,0);
% run_2state_blueWhale_IBM(6,3,2,0);
% run_2state_blueWhale_IBM(6,3,1,0);
% run_2state_blueWhale_IBM(6,3,4,70);
% run_2state_blueWhale_IBM(6,3,2,70);
% run_2state_blueWhale_IBM(6,3,1,70);
% run_2state_blueWhale_IBM(9,3,4,0);
% run_2state_blueWhale_IBM(9,3,2,0);
% run_2state_blueWhale_IBM(9,3,1,0);
% run_2state_blueWhale_IBM(9,3,4,70);
% run_2state_blueWhale_IBM(9,3,2,70);
% run_2state_blueWhale_IBM(9,3,1,70);
% run_2state_blueWhale_IBM(12,3,4,0);
% run_2state_blueWhale_IBM(12,3,2,0);
% run_2state_blueWhale_IBM(12,3,1,0);
% run_2state_blueWhale_IBM(12,3,4,70);
% run_2state_blueWhale_IBM(12,3,2,70);
% run_2state_blueWhale_IBM(12,3,1,70);

% run_2state_blueWhale_IBM(3,8,4,0);
% run_2state_blueWhale_IBM(3,8,2,0);
% run_2state_blueWhale_IBM(3,8,1,0);
% run_2state_blueWhale_IBM(3,8,4,70);
% run_2state_blueWhale_IBM(3,8,2,70);
% run_2state_blueWhale_IBM(3,8,1,70);
% run_2state_blueWhale_IBM(6,8,4,0);
% run_2state_blueWhale_IBM(6,8,2,0);
% run_2state_blueWhale_IBM(6,8,1,0);
% run_2state_blueWhale_IBM(6,8,4,70);
% run_2state_blueWhale_IBM(6,8,2,70);
% run_2state_blueWhale_IBM(6,8,1,70);
% run_2state_blueWhale_IBM(9,8,4,0);
% run_2state_blueWhale_IBM(9,8,2,0);
% run_2state_blueWhale_IBM(9,8,1,0);
% run_2state_blueWhale_IBM(9,8,4,70);
% run_2state_blueWhale_IBM(9,8,2,70);
% run_2state_blueWhale_IBM(9,8,1,70);
% run_2state_blueWhale_IBM(12,8,4,0);
% run_2state_blueWhale_IBM(12,8,2,0);
% run_2state_blueWhale_IBM(12,8,1,0);
% run_2state_blueWhale_IBM(12,8,4,70);
% run_2state_blueWhale_IBM(12,8,2,70);
% run_2state_blueWhale_IBM(12,8,1,70);

% create clouds
% create_clouds(3,1)
% create_clouds(3,3)
% create_clouds(3,8)
% create_clouds(6,1)
% create_clouds(6,3)
% create_clouds(6,8)
% create_clouds(9,1)
% create_clouds(9,3)
% create_clouds(9,8)
% create_clouds(12,1)
% create_clouds(12,3)
% create_clouds(12,8)

% close all;