function alpha = alpha(lambda, muest, M, N)
    alpha = zeros(M, 1);
    Ik = Ik_f(lambda, muest, M, N);
    for k = 1:M
        alpha(k) = alpha_k(lambda, muest, M, N, k) - 2*Ik;
    end

end

% function alpha_k = alpha_k(lambda, muest, M, N, k)
%     aux1 = 0;
%     aux2 = 0;
% 
%     aux3 = 0;
%     aux4 = 0;
%     aux5 = 0;
%     for r = 1:M
%         if r ~= k
%             % aux 1 
%             aux1 = aux1 + ( lambda(r) / (lambda(r) - lambda(k)) )*log(lambda(r) - 1)^2;
%             if aux1 < 0
%                  o = 1;
%             end
%             % aux 3
%             aux3 = aux3 + lambda(r) / (lambda(r) - lambda(k)); 
% 
%             aux4 = aux4 + log(lambda(k) / lambda(r))*log(lambda(k)/ abs(lambda(r) - lambda(k)));
%         end
%         % aux 2
%         aux2 = aux2 + ( muest(r) / (muest(r) - lambda(k))  ) * (log(muest(r)) - 1)^2;
%         aux3 = aux3 - lambda(r) / (lambda(r) - muest(k));
%         aux4 = aux4 - log(lambda(k) / muest(r))*log(lambda(k) / abs(muest(r) - lambda(k)));
% 
%         aux5 = aux5 + dilog(1 - (lambda(r) / lambda(k))) - dilog(1 - (muest(r) / lambda(k)));
%     end
% 
%     alpha_k = aux1 - aux2 + aux3*(log(lambda(k)) - 1)^2 + 2*log(lambda(k))^2 - 1 -2*aux4 - 2*aux5;
% 
%     if alpha_k < 0
%         display("here")
%     end
% end

function alpha_k = alpha_k(lambda, muest, M, N, k)
    lambda_less = lambda([1:k-1 k+1:M]);

    alpha_k_other1 = sum( (log(muest).^2) .* (muest ./ (lambda(k) - muest))  );
    alpha_k_other2 = -sum( (log(lambda_less).^2) .* (lambda_less ./ (lambda(k) - lambda_less))    );

    alpha_k_other3 = + 2*log(lambda(k)).^2 + 2*log(lambda(k)); %- 2*Ik_f(lambda, muest, M, N);

    alpha_k_other4 = - (log(lambda(k))^2)*( sum( lambda ./ (lambda - muest(k)) ) - sum(lambda_less ./ (lambda_less - lambda(k))) );


    alpha_k = alpha_k_other1 + alpha_k_other2 + alpha_k_other3 + alpha_k_other4;

%     alpha_k1 = 0;
%     alpha_k2 = 0;
%     alpha_k3 = 0;
%     alpha_k4 = 0;

%     alpha_k = 0;
% 
%     for r = 1:M
%         alpha_k1 = alpha_k1 + ((log(muest(r)))^2) * muest(r)/(lambda(k) - muest(r)) ;
%         if r ~= k
%             alpha_k2 = alpha_k2  - ((log(lambda(r)))^2) *lambda(r) / (lambda(k) - lambda(r));
%             alpha_k4 = alpha_k4 + (log(lambda(k))^2)*(lambda(r) / (lambda(r) - lambda(k)));
%         end
%         alpha_k4 = alpha_k4 - (log(lambda(k))^2)*(lambda(r) / (lambda(r) - muest(k)));
%     end
%     alpha_k3 =  2*log(lambda(k))^2 + 2*log(lambda(k)) - 2*Ik_f(lambda, muest, k, M, N);
% 
%     alpha_k = alpha_k1 + alpha_k2 + alpha_k3 + alpha_k4 ;
end



% 
% function alpha = alpha(lambda, muest, M,N)
% 
% aux1 = ((N-M)/M).*( ( log(muest) - 1 ).^2 - (log(lambda) - 1).^2 );
% aux2 = (1/M)*sum( log(lambda).^2 + 2*log(lambda) ) - 2;
% 
% 
% aux3 = 0;
% aux4 = 0;
% for k = 1:M
%     for r = 1:M
%         if r ~= k
%             aux3 = aux3 + log(lambda(k)  / lambda(r))*log( lambda(k)/abs(lambda(r) - lambda(k)) );
%         end
%         aux3 = aux3 - log(lambda(k) / muest(r))*log(lambda(k) / abs(muest(r) - lambda(k)));
% 
%         aux4 = aux4 + dilog(1 - (lambda(r) / lambda(k))) - dilog(1 - (muest(r)/ lambda(k)));
%     end
% end
% 
% alpha = sum(aux1) + aux2 - aux3 - real(aux4);
% 
% 
% end