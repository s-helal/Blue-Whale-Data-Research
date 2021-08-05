close all;

% gamma functions
mu = [3700, 1050];
sigma = [2150, 970];
step_length = 0:5:125000; % step lengths in meters
sst = 6:0.1:25; % degrees celcius
% for transit state 0:5:125000
% for forage state 0:5:100000

figure(1); 
figure(2); 
for rate=[1,2,4]
    [prob_t,prob_f] = get_dist(step_length,mu,sigma,rate);
    figure(1);
    plot(step_length, prob_t); hold on;
    figure(2);
    plot(step_length, prob_f); hold on;
end

figure(1); 
% title("Step Length Distribution in Forage State");
legend('Rate 1', 'Rate 2', 'Rate 4');
ylabel("Probability");
xlabel("Step Length");

xlim([0 125000]);
hold off;
saveas(gcf, "transit_pdf.png");
figure(2); 
% title("Turning Angle Distribition in Forage State");
legend('Rate 1', 'Rate 2', 'Rate 4');
xlim([0 125000]);
ylabel("Probability");
xlabel("Step Length");
hold off;
saveas(gcf, "forage_pdf.png");

% logit functions
% givens
alpha_1 = -1;
alpha_2 = -0.2;
x_star = 16;
beta = 5;
rho_star = 0.3;
% dists and data
x = sst;
rho = 0:0.1:1;
invlogit = @(x) 1./(1+exp(-x));
p_sst = invlogit(alpha_1 + alpha_2*(x - x_star));
p_krill = invlogit(beta*(rho - rho_star));
% plot
figure(3);
plot(x, p_sst);
ylabel("Probability of Whale in Forage State Based on SST")
xlabel("SST")
xlim([6 25]);
saveas(gcf, "change_sst_pdf.png");
figure(4);
plot(rho, p_krill);
ylabel("Probability of Whale in Forage State Based on Krill")
xlabel("Krill Density")
xlim([0 1]);
saveas(gcf, "change_krill_pdf.png");
