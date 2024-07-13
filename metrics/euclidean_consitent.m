function out = euclidean_consitent()
    
    out.dist = @euclidean_consitent_dist;
    out.plugin_dist = @euclidean_plug_in_dist;

    out.second_mean = @euclidean_consitent_second_mean;
    out.fluctuation_var = @euclidean_consistnet_fluctuation_var;

end

function dist = euclidean_plug_in_dist(r1, r2)
    M = size(r1, 1);
    
    dist = trace( (r1 -r2 )^2 )/M;
end

function dist = euclidean_consitent_dist(r1, r2, n1, n2)
    M = size(r1, 1);
    
    dist = trace( (r1 -r2 )^2 );
    
    dist = dist - (1/(n1))*trace(r1)^2 - (1/(n2))*trace(r2)^2;

    dist = dist / M;
end


function mean = euclidean_consitent_second_mean(r1, r2, n1, n2)
    
    M = size(r1, 1);
    
    mean = (1/(n1))*trace(r1^2) - (1/(n2))*trace(r2^2);
    
    % mean = mean / M;
end


function fluct_variance = euclidean_consistnet_fluctuation_var(r1, r2, n1, n2)
     M = size(r1, 1);
        
     Delta = r1 - r2;

     fluct_variance = 2*( (1/n1)*trace(r1^2) )^2 + 4*(1/n1)*trace( r1*Delta*r1*Delta );
     fluct_variance = fluct_variance + 2*( (1/n2)*trace(r2^2) )^2 + 4*(1/n2)*trace( r2*Delta*r2*Delta );
     fluct_variance = fluct_variance + 4*(1/(n1*n2))*trace(r1*r2)^2;

    
     % fluct_variance = fluct_variance/(M);


    fluct_variance = fluct_variance*2;  % real-valued
end