function out = fisher_consitent()
    % Code adapted from https://github.com/maliktiomoko/RMTEstimCovDist
    out.dist = @fisher_consitent_dist;
    out.plugin_dist = @fisher_plug_in_dist;

    out.second_mean = 0;
    out.fluctuation_var = 0;

end

function dist_trad = fisher_plug_in_dist(hatC1, hatC2)
    F=hatC1\hatC2;
    lambda=sort(eig(F));
    f = @(t) log(t).^2;
    dist_trad = mean(f(lambda));
end

function dist_consistent = fisher_consitent_dist(hatC1, hatC2, n1, n2)
    p = size(hatC1, 1);
    
    c2=p/n2;
    c1=p/n1;                   

    F=hatC1\hatC2;
    
    lambda=sort(eig(F));
    slambda=sqrt(lambda);
    eta = sort(eig(diag(lambda)-slambda*slambda'/(p-n1)));
    zeta = sort(eig(diag(lambda)-slambda*slambda'/n2));
    
    
    
    M=zeros(p);
    N=zeros(p);                
    for i=1:p
        M(i,i) = 1/(2*lambda(i)^2);
        N(i,i) = 1/lambda(i);                 
        js = 1:p;
        js(i) = [];
        for j=js
            M(i,j) = (-1+lambda(i)/lambda(j)-log(lambda(i)/lambda(j)))/(lambda(i)-lambda(j))^2;                        
            N(i,j) = log(lambda(i)/lambda(j))/(lambda(i)-lambda(j));
        end
    end     
    
    %%% Large p-estimate
    
    dist_consistent = 2*(c1+c2-c1*c2)/(c1*c2)*( (eta-zeta)'*M*(eta-lambda)+(eta-lambda)'*(log((1-c1)*lambda)./lambda) )...
                     -2/p*(eta-zeta)'*N*ones(p,1)+1/p*sum(log((1-c1)*lambda).^2)...
                     -2*(1-c2)/c2*( 1/2*log( (1-c1)*(1-c2) )^2+(eta-zeta)'*(log((1-c1)*lambda)./lambda) );  
    
    
    
end



