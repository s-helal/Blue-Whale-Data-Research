function [t_axis,mean_vals,std_vals_minus,std_vals_plus] = make_mean_plots(percent)%(n,m,steps_per_day)
close all;
%     filename1 = 'input_data_geo_' + string(n) + '_tmp_' + string(m) + '_yr_2008.mat';
%     filename2 = ['output_data_geo_' num2str(n) '_tmp_' num2str(m) '_rate_' num2str(steps_per_day) '_yr_2008.mat'];  
    filename1 = ['not_holey_input_' num2str(percent) '_percent.mat'];
    filename2 = ['not_holey_output_' num2str(percent) '_percent.mat'];  

    load(filename1);
    load(filename2);
    
    
    for k = 1:size(Y,1)*size(Y,2)
        if isnan(Y(k))
        Y(k) = 0;
        end
    end
%      
%     nc_file_names = dir('Data/wc15n_srf_2008.nc');
% 
%     % From run_2state_blueWhale_IBM.m
%     nc_file = [nc_file_names.folder '/' nc_file_names.name];
%     
%     % Code snip-it to turn IBM (X,Y) coordinates (in meters) to
%     % latitude-longitude coordinates.
%     longit = ncread(nc_file, 'lon_rho');  longit = permute(longit,[2,1]);
%     latit  = ncread(nc_file, 'lat_rho');  latit  = permute(latit, [2,1]);
% 
%     %long_ref = min(longit(:)); % our (0,0)
%     lat_ref  = min(latit(:));
% 
%     rad_earth = 6378.137;  % Radius of the earth (KM)
%     m1 = 360/(2*pi*rad_earth*1000); % Scaling factor
% 
%     %Map coordinates to lat-long
%     Y_lat  = lat_ref + Y.*m1;
%     %X_long = long_ref + (X.*m1)./cos(Y_lat.*pi./180);
%  
%     % Monterey Bay Latitudes
%     MB_lat_S = 36.6; % Southern boundary
%     MB_lat_N = 37;   % Northern boundary
%     
    % initialize y-axis vectors
      mean_vals = zeros(size(Y,2),1);
      std_vals_minus = zeros(size(Y,2),1);
      std_vals_plus = zeros(size(Y,2),1);

        %mode_vals = zeros(size(Y,2),1);
        state_vals = zeros(size(s,2),1);
%         monterey_vals = zeros(size(Y,2),1);
    
        q1_vals = zeros(size(Y,2),1);
        median_vals = zeros(size(Y,2),1);
        q3_vals = zeros(size(Y,2),1);

    for k = 1:size(Y,2)
        Y_new = nonzeros(Y(:,k)); 
        mean_vals(k) = mean(Y_new);
        std_Y = std(Y_new, 0, 'all', 'omitnan');
        std_vals_minus(k) = mean(Y_new)-std_Y;
        std_vals_plus(k) = mean(Y_new)+std_Y;
        %mode_vals(k) = mode(Y_new);
        state_vals(k) = size(nonzeros(s(:,k) == 2),1)/size(nonzeros(s(:,k)),1);
%         monterey_vals(k) = size(nonzeros((Y_lat(:,k) >= MB_lat_S) & (Y_lat(:,k) <= MB_lat_N)),1);

        q1_vals(k) = median(Y_new(Y_new > median(Y_new)));
        median_vals(k) = median(Y_new);
        q3_vals(k) = median(Y_new(Y_new < median(Y_new)));
    end

    % t-axis is one unit = one timestep
     t_axis = linspace(1,par.numDays,par.numDays*par.rate); %1:par.numDays*par.rate;
%     title_tail = [num2str(n) ' km, ' num2str(m) ' days, ' num2str(steps_per_day)  ' steps/day'];
    title_tail = [num2str(percent) '%'];

    % figures
    figure(1);
    plot(t_axis,mean_vals); xlabel('timestep'); ylabel('mean Y');
    title(['Mean Latitude (' title_tail ')'] ); set(gca,'fontsize',10);

    %figure(2);
    %plot(t_axis,mode_vals); xlabel('timestep'); ylabel('mode Y');
    %title('Mode Latitudinal Distribution vs. Time');

    figure(3);
    plot(t_axis,state_vals);
    title(['Travel State (' title_tail ')'] ); set(gca,'fontsize',10);

%     figure(4);
%     plot(t_axis,monterey_vals); xlabel('timestep'); ylabel('whales in Monterey Bay');
%     title(['Whales in Monterey (' title_tail ')'] ); set(gca,'fontsize',10);

    figure(5); hold on;
    xlabel('timestep'); ylabel('latitude');
    title(['Quantile Latitudes (' title_tail ')'] ); set(gca,'fontsize',10);
%     p1 = 
    plot(t_axis,q1_vals); %a1 = "Q1";
%     p2 = 
    plot(t_axis,median_vals); %a2 = "median";
%     p3 = 
    plot(t_axis,q3_vals); %a3 = "Q3";
    %legend([p1,p2,p3],[a1,a2,a3]);
    hold off

%     name_tail = [num2str(n) '_tmp_' num2str(m) '_yr_2008_rate_' num2str(steps_per_day) '.png'];
%     saveas(figure(1), ['mean_plot_geo_' name_tail])
%     %saveas(figure(2), 'lat_mode_plot_geo_' + name_tail)
%     saveas(figure(3), ['state_plot_geo_' name_tail])
%     saveas(figure(4), 'monterey_plot_geo_' + name_tail)
%     saveas(figure(5), ['quantile_plot_geo_' name_tail])

    saveas(figure(1), ['not_holey_mean_' num2str(percent) '_percent.png'])
    saveas(figure(3), ['not_holey_state_' num2str(percent) '_percent.png'])
    saveas(figure(5), ['not_holey_quantile_' num2str(percent) '_percent.png'])

% figure
% [ta,a,am,ap] = make_mean_plots(3,1,1);
% [tb,b,bm,bp] = make_mean_plots(3,1,2);
% [tc,c,cm,cp] = make_mean_plots(3,1,4);
% plot(ta,a,'r',ta,am,'-r',ta,ap,'-r',tb,b,'k',tb,bm,'-k',tb,bp,'-k',tc,c,'b',tc,cm,'-b',tc,cp,'-b');
% title('3 km, 1 day, 1-2-4 steps: Mean and Stdv')

%% get percent foraging for clouded-- run alone
    % with clouds n=3,m=1; n=3,m=3; n=12,m=1; n=12,m=3
    % without clouds all above + n=3,m=8; n=12; m=8
    n=3;m=1;
    inpath = 'C:\Users\samih\OneDrive\Documents\Blue Whale Data Resolution Project\New Clouded Data\';
    name_id_cloud = ['geo_' num2str(n) '_tmp_' num2str(m) '_clouded_' num2str(8)];
    cloud_file = [inpath 'output_data_' name_id_cloud '.mat'];
    load(cloud_file)
    
    for k = 1:size(Y,1)*size(Y,2)
        if isnan(Y(k))
        Y(k) = 0;
        end
    end
    state = zeros(size(s,2),1);
    for k=1:size(Y,2)
        state(k) = size(nonzeros(s(:,k) == 2),1)/size(nonzeros(s(:,k)),1);
    end
    writematrix(state, [inpath 'state_' name_id_cloud '.csv']);
    
    inpath2 = 'C:\Users\samih\OneDrive\Documents\Blue Whale Data Resolution Project\New Output Data\';
    name_id_reg = ['geo_' num2str(n) '_tmp_' num2str(m) '_yr_2008'];
    reg_file = [inpath2 'output_data_' name_id_reg '.mat'];
    load(reg_file)
    
    for k = 1:size(Y,1)*size(Y,2)
        if isnan(Y(k))
        Y(k) = 0;
        end
    end
    state = zeros(size(s,2),1);
    for k=1:size(Y,2)
        state(k) = size(nonzeros(s(:,k) == 2),1)/size(nonzeros(s(:,k)),1);
    end
    writematrix(state, [inpath 'state_' name_id_reg '.csv']);
end