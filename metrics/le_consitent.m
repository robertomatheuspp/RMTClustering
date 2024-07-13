function out = le_consitent()
    
    out.dist = @le_consitent_dist;
    out.plugin_dist = @le_plug_in_dist;

    out.second_mean = @le_consitent_second_mean;
    out.fluctuation_var = @le_consistnet_fluctuation_var;

end

function dist = le_plug_in_dist(r1, r2)
    M = size(r1, 1);
    
    dist = trace( (logm(r1) - logm(r2))^2 ) / M;
    
end

function dist = le_consitent_dist(r1, r2, n1, n2)
    M = size(r1, 1);
    
    [dist, ~, ~] = le_imp(r1, r2, n1, n2);    
end


function mean = le_consitent_second_mean(r1, r2, n1, n2)
    % M = size(r1, 1);
    % 
    % [eigv1, lamb1] = eig(r1); lamb1 = diag(lamb1);
    % [lamb1, idx_sort1] = sort(lamb1); 
    % eigv1 = eigv1(:, idx_sort1);
    % 
    % 
    % [eigv2, lamb2] = eig(r2); lamb2 = diag(lamb2);
    % [lamb2, idx_sort2] = sort(lamb2);
    % eigv2 = eigv2(:, idx_sort2);
    % 
    % 
    % [theta1_aux, theta1] = rooting_theta(lamb1, M, n1);
    % [theta2_aux, theta2] = rooting_theta(lamb2, M, n2);
    % 
    % aux1 = 0; aux2 = 0; aux3 = 0; aux4 = 0;
    % 
    % for m = 1:M
    %     proj_m1 = eigv1(m, :)'*eigv1(m, :);
    %     proj_m2 = eigv2(m, :)'*eigv2(m, :);
    % 
    %     aux1 = aux1 - trace(proj_m1*(logm(r2)  - log(lamb1(m))*eye(M) )^2);
    %     aux2 = aux2 - trace(proj_m2*(logm(r1)  - log(lamb2(m))*eye(M) )^2);
    % 
    % 
    % 
    %     aux3 
    % end
    % 

end


function fluct_variance = le_consistnet_fluctuation_var(r1, r2, n1, n2)
    M = size(r1, 1);
    
    fluct_variance =  leVariance(r1, r2, M, n1, n2);
end


