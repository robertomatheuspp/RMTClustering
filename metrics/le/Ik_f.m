% function Ik = Ik_f(lambda, muest, k, M, N)
% 
%     lambda_less = lambda([1:k-1 k+1:M]);
% 
%     Ik = 1 - log(N)*log(lambda(k)) + log(lambda(k))^2 + (log(lambda(k)) - 1)*log(1 - M/N);
%     Ik = Ik - sum( (lambda(k)./ (lambda(k) - lambda_less)).*log(lambda_less)  );
% 
%     Ik = Ik + sum( (lambda(k)./ (lambda(k) - muest) ) .*log(muest) );
% 
%     Ik = Ik + log(lambda(k))*( sum(  lambda_less ./ (lambda_less - lambda(k)) )  -   sum(lambda./(lambda - muest(k))) );
% 
%     Ik = Ik + sum( log(lambda_less).*log(abs(lambda(k) - lambda_less )) ) - sum( log(muest).*log(abs(lambda(k) - muest)) );
%     Ik = Ik + sum(dilog(1 - lambda./lambda(k)) - dilog( 1 - muest ./ lambda(k)));
%     
%     Ik = real(Ik);



function Ik = Ik_f(lambda, muest, M, N)

    Ik = 0.5*(N/M - 1)*log(1-M/N)^2 - (N/M - 1)*log(1 - M/N) - 1;

    aux1 = 0;
    aux2 = 0;
    for k = 1:M
        lambda_less = lambda([1:k-1 k+1:M]);
        aux1 = aux1 + sum(  log( muest ).*log( lambda(k) ./ (lambda(k) - muest) )  );
        aux1 = aux1 - sum(  log(lambda_less).*log( lambda(k) ./ (lambda(k) - lambda_less) )  );

        aux2 = aux2 + sum(  phi2(lambda./lambda(k)) - phi2(muest./lambda(k))  );
    end

    aux1 = aux1 / M;
    aux2 = aux2 / M;

    Ik = Ik + aux1 + aux2;

    Ik = Ik - log(N)*(1/M)*sum(  log(lambda)  );

    Ik = real(Ik);
end




function out = phi2(x)  
    out = zeros(length(x), 1);
    out(x < 1)  =  dilog(1 - x(x < 1));
    out(x >= 1) = (pi^2)/3 - (1/2)*(log(x(x >= 1)).^2) - dilog(1 - 1./x(x >= 1));
    %     if x < 1
    %         out = dilog(1 - x);
    %     else
    %         out = (pi^2)/3 - (1/2)*(log(x).^2) - dilog(1 - 1./x);
    %     end
    % 
end