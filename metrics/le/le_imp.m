function [dist_le, alpha1_mean, alpha2_mean] = le_imp(cov1est, cov2est, N1, N2)
    M = size(cov1est, 1);

    [E1, lambda1] = eig(cov1est);
    [lambda1, idx1] = sort(diag(lambda1));
    E1 = E1(:, idx1);
    muest1 = newton_rapson(lambda1, M/N1);

    [E2, lambda2] = eig(cov2est);
    [lambda2, idx2] = sort(diag(lambda2));
    E2 = E2(:, idx2);
    muest2 = newton_rapson(lambda2, M/N2);  


    beta1 = beta(lambda1, muest1);
    beta2 = beta(lambda2, muest2);

    alpha1_mean = mean(alpha2(lambda1, muest1, M, N1));
    alpha2_mean = mean(alpha2(lambda2, muest2, M, N2));


    dist_le = alpha1_mean + alpha2_mean - (2/M)*trace( (E1*diag(beta1)*E1') * (E2*diag(beta2)*E2') );
    

    if isnan(dist_le) 
        display("here")
    end

    if isinf(dist_le) 
        display("here")
    end
    
end