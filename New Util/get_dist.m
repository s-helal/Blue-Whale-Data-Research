function [prob_t, prob_f] = get_dist(data, mu, sigma, rate)
    % givens
    hours = 24/rate;
    mu = hours*mu; % [transit, forage]
    sigma = hours*sigma;
    % get parameters 
    alpha_t = mu(1)^2 / sigma(1)^2;
    beta_t = mu(1) / sigma(1)^2;
    alpha_f = mu(2)^2 / sigma(2)^2;
    beta_f = mu(2) / sigma(2)^2;
    % get distribution and data
    g = @(x,alpha,beta) beta^alpha*x.^(alpha-1).*exp(-beta .* x)./gamma(alpha);
    prob_t = g(data,alpha_t,beta_t);
    prob_f = g(data,alpha_f,beta_f);
end
