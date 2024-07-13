function betak = beta(lambda, muest)
M = size(lambda, 1);
% betak=zeros(M, 1); % Weights of the proposed method
% for m=1:M,
%     betak(m)=1;
%     temp1=0;
%     for k=1:M,
%         if k~=m,
%             temp1=temp1+lambda(m)/(lambda(k)-lambda(m));
%         end;
%     end;
%     betak(m)=betak(m)+(1+temp1-sum(muest(m)./(lambda-muest(m))))*log(lambda(m));
%     temp2=0;
%     for r=1:M,
%         if r~=m,
%             temp2=temp2+lambda(r)/(lambda(r)-lambda(m))*log(lambda(r));
%         end;
%     end;
%     betak(m)=betak(m)+temp2+sum(muest./(lambda(m)-muest).*log(muest));
% end;



betak = zeros(M, 1);
for k = 1:M
    lambda_less = lambda([1:k-1 k+1:M]);

    betak(k) = 1 + (1 + sum( lambda(k) ./ (lambda_less - lambda(k)) ) - sum( muest(k) ./ (lambda - muest(k)) ) ).*log(lambda(k));
    betak(k) = betak(k) + sum(log(lambda_less) .* (lambda_less ./(lambda_less - lambda(k)) ) ) - sum(( muest./(muest - lambda(k))).*log(muest) );

end

end