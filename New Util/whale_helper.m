function whale_helper(percent)% percent = 30;    
orig_file = 'output_data_geo_3_tmp_1_rate_4_yr_2008.mat';
    edit_file = ['not_holey_output_' num2str(percent) '_percent.mat'];
    load(orig_file);
    Y_orig = Y;
    load(edit_file, 'Y');
    Y_edit = Y;
    
    for k = 1:size(Y_orig,1)*size(Y_orig,2)
        if isnan(Y_orig(k))
        Y_orig(k) = 0;
        end
    end
    
    for k = 1:size(Y_edit,1)*size(Y_edit,2)
        if isnan(Y_edit(k))
        Y_edit(k) = 0;
        end
    end
    
    Y_diff = zeros(size(Y,2),1);
    for k = 1:size(Y,2)
        Y_orig_new = nonzeros(Y_orig(:,k));
        Y_edit_new = nonzeros(Y_edit(:,k));
        Y_diff(k) = mean(Y_orig_new)-mean(Y_edit_new);

    end
    
    std_Y = std(Y_diff, 0, 'all', 'omitnan');
    std_vals_minus = Y_diff-std_Y;
    std_vals_plus = Y_diff+std_Y;
    
    t = linspace(1,par.numDays,par.numDays*par.rate); %1:par.numDays*par.rate;
    plot(t, Y_diff, t, std_vals_minus, '-', t, std_vals_plus,'-');
    title(['Mean Location Difference ' num2str(percent) ' Percent']);
    saveas(gcf, ['not_holey_diff_' num2str(percent) '_percent.png']);
end
    
% GENERATE COMMAND for run_2state_blueWhale_IBM

% n_array = [3,9,24];
% m_array = [1,8];
% steps_array = [1,2,4];
% 
% for i = 1:3
%     for j = 1:2
%         for k = 1:3
%             name = ['n = ' num2str(n_array(i)) '; m = ' num2str(m_array(j)) '; steps_per_day = ' num2str(steps_array(k)) '; run_2state_blueWhale_IBM; clear all;'];
%             disp(name)
%         end
%     end
% end

% n = 3; m = 1; steps_per_day = 1; run_2state_blueWhale_IBM; clear all;
% n = 3; m = 1; steps_per_day = 2; run_2state_blueWhale_IBM; clear all;
% n = 3; m = 1; steps_per_day = 4; run_2state_blueWhale_IBM; clear all;
% n = 3; m = 8; steps_per_day = 1; run_2state_blueWhale_IBM; clear all;
% n = 3; m = 8; steps_per_day = 2; run_2state_blueWhale_IBM; clear all;
% n = 3; m = 8; steps_per_day = 4; run_2state_blueWhale_IBM; clear all;
% n = 9; m = 1; steps_per_day = 1; run_2state_blueWhale_IBM; clear all;
% n = 9; m = 1; steps_per_day = 2; run_2state_blueWhale_IBM; clear all;
% n = 9; m = 1; steps_per_day = 4; run_2state_blueWhale_IBM; clear all;
% n = 9; m = 8; steps_per_day = 1; run_2state_blueWhale_IBM; clear all;
% n = 9; m = 8; steps_per_day = 2; run_2state_blueWhale_IBM; clear all;
% n = 9; m = 8; steps_per_day = 4; run_2state_blueWhale_IBM; clear all;
% n = 24; m = 1; steps_per_day = 1; run_2state_blueWhale_IBM; clear all;
% n = 24; m = 1; steps_per_day = 2; run_2state_blueWhale_IBM; clear all;
% n = 24; m = 1; steps_per_day = 4; run_2state_blueWhale_IBM; clear all;
% n = 24; m = 8; steps_per_day = 1; run_2state_blueWhale_IBM; clear all;
% n = 24; m = 8; steps_per_day = 2; run_2state_blueWhale_IBM; clear all;
% n = 24; m = 8; steps_per_day = 4; run_2state_blueWhale_IBM; clear all;

% GENERATE COMMAND for make_mean_plots
% for i = 1:3
%     for j = 1:2
%         for k = 1:3
%             name = ['make_mean_plots(' num2str(n_array(i)) ',' num2str(m_array(j)) ',' num2str(steps_array(k)) '); close all; clear all;'];
%             disp(name)
%         end
%     end
% end

% make_mean_plots(3,1,1); close all; clear all;
% make_mean_plots(3,1,2); close all; clear all;
% make_mean_plots(3,1,4); close all; clear all;
% make_mean_plots(3,8,1); close all; clear all;
% make_mean_plots(3,8,2); close all; clear all;
% make_mean_plots(3,8,4); close all; clear all;
% make_mean_plots(9,1,1); close all; clear all;
% make_mean_plots(9,1,2); close all; clear all;
% make_mean_plots(9,1,4); close all; clear all;
% make_mean_plots(9,8,1); close all; clear all;
% make_mean_plots(9,8,2); close all; clear all;
% make_mean_plots(9,8,4); close all; clear all;
% make_mean_plots(24,1,1); close all; clear all;
% make_mean_plots(24,1,2); close all; clear all;
% make_mean_plots(24,1,4); close all; clear all;
% make_mean_plots(24,8,1); close all; clear all;
% make_mean_plots(24,8,2); close all; clear all;
% make_mean_plots(24,8,4); close all; clear all;
