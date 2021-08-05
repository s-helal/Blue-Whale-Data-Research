filename = 'erd_GOES_sst_1day_2019_formatted.mat';
load(filename);
cover = zeros(size(sst,3),1);
for i = 1:size(sst,3)
    cover(i) = numel(nonzeros(sst(:,:,i)<0));
end
ratio = cover/(nnz(~isnan(sst(:,:,1))));
figure;
hold on;
plot(1:size(sst,3),ratio,'LineWidth',1.25);
title('GOES 1 Day Percent Cloud Cover');
for j=0:31:364
    plot(j*ones(100),linspace(-0.05,0.35),'--')
    % xline(j)
end
yline(0);
hold off;
saveas(gcf, 'goes_1d_cloudpercent.png');