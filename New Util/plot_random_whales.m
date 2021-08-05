function plot_random_whales(n,m)
    % plot some random whales: copied from the IBM driver
    filename1 = 'input_data_geo_' + string(n) + '_tmp_' + string(m) + '_yr_2008.mat';
    filename2 = 'output_data_geo_' + string(n) + '_tmp_' + string(m) + '_yr_2008.mat';
    load(filename1);
    load(filename2);

    t_vals = 1:par.numDays*par.rate;
    cw = randi(size(Y,1),1,3);    % cw-- chosen whales
    whale_lats = Y(cw,:);

    % from main
    X_domain = grid_pars.xrange(1):grid_pars.resolution:grid_pars.xrange(2); X_domain = X_domain(1:end-1);
    Y_domain = grid_pars.yrange(1):grid_pars.resolution:grid_pars.yrange(2); Y_domain = Y_domain(1:end-1);

    [XX, YY] = meshgrid(X_domain, Y_domain);

    % if all plots should stay open: figure(n);
    hold on;
    pcolor(XX,YY,krill_data(:,:,end)); shading interp;
    plot(X(cw,:)',Y(cw,:)', 'LineWidth', 3);
    hold off;
    saveas(gcf, char('rand_lat_plot_geo_' + string(n) + '_tmp_' + string(m) + '_yr_2008.png'))

    figure(2)
    plot(t_vals,whale_lats)
    saveas(gcf, char('rand_whale_map_geo_' + string(n) + '_tmp_' + string(m) + '_yr_2008.png'))
end