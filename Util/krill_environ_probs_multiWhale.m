function [gamma_K, gamma_E, gamma] = krill_environ_probs_multiWhale(krill_loc, sst_loc, good_whales, par)

% gamma_K: trans prob based on krill
% gamma_E: trans prob based on environment
% gamma  : trans prob based on combo of krill and environ
% Size of gamma matrices depends on krill_loc and sst_loc inputs 


% Rows are whales, columns are locations (multiple locations if conducting ARS)
N = size(sst_loc,2);
alpha = repmat(par.envir_pref(good_whales,1),1,N) + repmat(par.envir_pref(good_whales,2),1,N).*(sst_loc + repmat(par.envir_pref(good_whales,3),1,N));
gamma_E = exp(alpha)./(exp(alpha) + 1); 
gamma_K = (1 + exp(repmat(par.krill_pref(good_whales,1),1,N).*(-krill_loc + repmat(par.krill_pref(good_whales,2),1,N)))).^(-1);

gamma   = (par.w(1).* gamma_E + par.w(2).*gamma_K)./par.w(3);