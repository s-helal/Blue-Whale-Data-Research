function make_map(file_head)%(n,m,file_head)
    % where to get whale positions
    input_data_path = 'C:\Users\samih\OneDrive - University of California, Davis\Percent Input Data\Temp\';%Output Data\';
    % where to get krill and sst data
%     output_data_path = 'C:\Users\samih\OneDrive\Documents\Blue Whale Data Resolution Project\';
    % where to save images
    image_path = 'C:\Users\samih\OneDrive\Documents\Blue Whale Data Resolution Project\Misc Images\Temp\';
    
    file_head = replace(file_head, ["input_data_","output_data_",".mat"],"");
    % choose files
%     whale_file = ['output_data_' file_head '.mat'];
    map_file = ['input_data_' file_head '.mat'];
    
    
    % GET MAP
    load([input_data_path map_file], 'krill_data', 'sst_data');
%     X_domain = grid_pars.xrange(1):grid_pars.resolution:grid_pars.xrange(2); X_domain = X_domain(1:end-1);
%     Y_domain = grid_pars.yrange(1):grid_pars.resolution:grid_pars.yrange(2); Y_domain = Y_domain(1:end-1);

    % GET WHALES
%     load([output_data_path whale_file], 'X', 'Y', 'day_vec', 'grid_pars', 'par');

    % prepare to plot
%     [XX, YY] = meshgrid(X_domain, Y_domain);
    m = 1;
    if ~contains(file_head,"cloud") && contains(file_head,"tmp_8")
        m = 8;
    end
        
    timeslice = floor(245/m); % max(day_vec(:))
    % 245 means september, 8 is for 8 day, 1 for 1 day or clouded 1 day
    
    % save map and krill
%     writematrix([XX, YY], 'map.csv')
%     writematrix(krill_data(:,:,245), 'krill.csv')
    
    % PLOT MAP SST
%     figure; pcolor(XX,YY,sst_data(:,:,timeslice)); shading interp;  % if have whales/output data
    figure; pcolor(sst_data(:,:,timeslice)); shading interp;
    title(strrep(file_head, "_", " ") + " sst"); set(gca,'fontsize',14);
    saveas(gcf,[image_path 'map_' file_head '_sst' '.png']);
    
    % PLOT MAP KRILL
%     figure; pcolor(XX,YY,krill_data(:,:,timeslice)); shading interp;  % if have whales/output data
    figure; pcolor(krill_data(:,:,timeslice)); shading interp;

    title(strrep(file_head, "_", " ") + " krill"); set(gca,'fontsize',14);
    saveas(gcf,[image_path 'map_' file_head '_krill' '.png']);    
%     
%     % PLOT WHALES
%     hold on;
%     for j = 1:par.numWhales
%         plot(X(j,:),Y(j,:),'k.','linewidth',2);
%     end
% 
%     hold off; drawnow;
%     saveas(gcf,[outpath 'map_' file_head '_whales.png']);
%     
%     % print number of whales in domain
%     disp("number of whales in domain for " + file_head + " is " + sum(X(:,timeslice)>0))
   
end