function out = kl_consitent()
    
    out.dist = @kl_consitent_dist;
    out.plugin_dist = @kl_plug_in_dist;

    out.second_mean = @kl_consitent_second_mean;
    out.fluctuation_var = @kl_consistnet_fluctuation_var;

end

function dist = kl_plug_in_dist(r1, r2)
    M = size(r1, 1);
    
    dist = trace( (r2*pinv(r1) ))/(2*M)  +  trace( (r1*pinv(r2) ))/(2*M)  - 1;
    
end

function dist = kl_consitent_dist(r1, r2, n1, n2)
    M = size(r1, 1);
    
    dist = (1 - M/n1)*(1/(2*M))*trace( (pinv(r1)*r2 )) +  (1 - M/n2)*(1/(2*M))*trace( (pinv(r2)*(r1) ))  - 1;
    
    
end


function mean = kl_consitent_second_mean(r1, r2, n1, n2)
    M = size(r1, 1);

    mean = (1/(n1 - M))*trace(pinv(r1)*r2) + (1/(n2 - M))*trace(pinv(r2)*r1);

    mean = mean / 2;
end


function fluct_variance = kl_consistnet_fluctuation_var(r1, r2, n1, n2)
    M = size(r1, 1);
    
    pinv_r1_r2 = pinv(r1)*r2;
    r1_pinv_r2 = r1*pinv(r2);
    
    aux1 = -M/(2*n1*n2);
    
    aux2 = trace(r1_pinv_r2*r1_pinv_r2)/( 4*(n2 - M)*n1 );
    aux3 = trace(pinv_r1_r2*pinv_r1_r2)/( 4*(n1-M)*n2 );


    aux4 = (0.25/n1)*( trace(r1_pinv_r2)/(n2 - M) )^2;    
    aux5 = (0.25/n2)*( trace(pinv_r1_r2)/(n1 - M) )^2;    


    aux_mult =  2*(n1 + n2 - M); % 2 is for real-valued

    
    fluct_variance = aux_mult*(aux1 + aux2 + aux3 + aux4 + aux5);

end


