
nbsamples = 1e2;


indexSCM1 = [1 1 1 1 1 2 2 2 2 3 3 3 4 4 5];  % Contains i_1, i_2, ... , i_R
indexSCM2 = [2 3 4 5 6 3 4 5 6 4 5 6 5 6 6];  % Contains j_1, j_2, ... , j_R

% Index to consider when calculating INTRA cluster distance
idx_intra = [1, 10, 15];

% Index to consider when calculating INTER cluster distance
idx_inter = [2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14];




R = length(indexSCM1);



% scenario 1
clist = [1/2, 1/2, 1/3, 1/3, 1/4, 1/4];
rholist = [0.3, 0.3, 0.5, 0.5, 0.7, 0.7];

% scenario 2
% clist = [1/4, 1/4, 1/3, 1/3, 1/2, 1/2];
% rholist = [0.3, 0.3, 0.6, 0.6, 0.9, 0.9];


% Instanciating the functionals I want to use. Plugin and consistent are
% functions withing these instances.
metric_f_list = [euclidean_consitent(), ...
                 kl_consitent(), ...
                 le_consitent(), ...
                 fisher_consitent(), ...
                 wasserstein_consistent()];


% Storing distances so that later we can calulate the probability of
% success. dist(:,:,:, 1) will store plugin distance and 
% dist(:,:,:, 2) the consistent estimators 
dist = zeros(nbsamples, length(metric_f_list), R, 2);


Mrange = [10, 20, 30, 40];%, 50, 60, ];% 70, 80, 90, 100];

prob_succ = zeros(length(Mrange), length(metric_f_list), 2);

for idxM = 1:length(Mrange)
    M = Mrange(idxM);
    fprintf("Running M = %d\n", M)

    covMat_list = zeros(M, M, length(rholist)); 
    Nlist = round(M./clist);


    for idx_cov = 1:length(rholist)
        covMat_list(:, :, idx_cov) =  toeplitz(rholist(idx_cov).^(0:M-1)); 
    end
    for idx_sim = 1:nbsamples
        
        % generating samples
        curSCM = zeros(M, M, R);
        for r = 1:length(rholist)
            curY = sqrtm(covMat_list(:,:,r))*randn(M, Nlist(r));
            curSCM(:, :, r) = curY*curY' / Nlist(r);
        end
    
        % here we will store the result for this current Monte Carlo
        % exectution
        cur_dist = zeros(length(metric_f_list), R); 
        cur_pluging_dist = zeros(length(metric_f_list), R); 
         
        for idx_metric = 1:length(metric_f_list)
            cur_metric = metric_f_list(idx_metric);
            for r = 1:R
	            i_r = indexSCM1(r);
                j_r = indexSCM2(r);
                
            
                r1 = curSCM(:, :, i_r);
                r2 = curSCM(:, :, j_r);
                N1 = Nlist(i_r); 
                N2 = Nlist(j_r);
    
                cur_dist(idx_metric, r) =  cur_metric.dist(r1, r2, N1, N2);
                cur_pluging_dist(idx_metric, r) =  cur_metric.plugin_dist(r1, r2);
            end
        end
    
        % store distances
        dist(idx_sim, :, :, 1) = cur_pluging_dist; 
        dist(idx_sim, :, :, 2) = cur_dist; 
    end


    cur_prob_succ = max(dist(:, :, idx_intra, :), [], 3) < min(dist(:, :, idx_inter, :), [], 3);
    prob_succ(idxM, :, :) = mean(squeeze(cur_prob_succ), 1);

end

display("Monte Carlo finished.")


%%


figure(1)
hold on;
plot(Mrange, prob_succ(:, 1, 1), "r");
plot(Mrange, prob_succ(:, 1, 2), "r--");

plot(Mrange, prob_succ(:, 2, 1), "g");
plot(Mrange, prob_succ(:, 2, 2), "g--");

plot(Mrange, prob_succ(:, 3, 1), "b");
plot(Mrange, prob_succ(:, 3, 2), "b--");


plot(Mrange, prob_succ(:, 4, 1), "m");
plot(Mrange, prob_succ(:, 4, 2), "m--");

plot(Mrange, prob_succ(:, 5, 1), "k");
plot(Mrange, prob_succ(:, 5, 2), "k--");

hold off;


xlabel("M")
ylabel("Prob. Correct Clustering")
legend(["EU Plugin", "EU  Consistent", "KL Plugin", "KL  Consistent", "LE Plugin", "LE Consistent"]);


legend(["EU Plugin", "EU  Consistent", ...
        "KL Plugin", "KL  Consistent", ...
        "LE Plugin", "LE Consistent", ...
        "Fisher Plugin", "Fisher Consistent", ...
        "Wasserstein Plugin", "Wasserstein Consistent", ...
        ]);







