function alpha = alpha2(lambda, muest, M, N)

    aux1 = (N/M - 1).*log(muest).^2 - (N/M - 1).*log(lambda).^2;

    aux2 = sum( (1/M)*log(lambda).^2 ) + sum( (2/M).*log(lambda) ) + 2;
    aux3 = -(N/M -1)*log(1 - M/N)^2 + 2*(N/M - 1)*log(1 - M/N);

    aux4_1 = log(lambda' ./ lambda).*(log(lambda ./ abs(lambda'- lambda)) );
    aux4_1(isnan(aux4_1)) = 0;
    aux4_1 = sum(sum(aux4_1));

    aux4_2 = -sum( sum( log(muest ./ lambda' ).*log(lambda' ./ abs(muest - lambda' )) ) );
    aux5 = sum(phi2(muest./lambda')  - phi2(lambda./lambda'));

    alpha  = sum(aux1) + sum(aux2)  + aux3 + 2*(aux4_1 + aux4_2 + aux5)/M;
end



function out = phi2(x)  
    x = real(x);
    out = zeros(length(x), 1);

    % out(x < 1)  =  dilog(1 - x(x < 1));
    % out(x >= 1) = (pi^2)/3 - (1/2)*(log(x(x >= 1)).^2) - dilog(1 - 1./x(x >= 1)); 

    out(x < 1)  =  myLi2(x(x < 1));
    out(x >= 1) = (pi^2)/3 - (1/2)*(log(x(x >= 1)).^2) - myLi2(1./x(x >= 1)); 


end


function z = myLi2(x)
    fun = @(x) log(1 - x) ./ x;
    
    
    % trapz(X,Y)
    aux = linspace(0.0000001, 0.999999, 500);
    xaux = x*aux;
    delta = xaux(:, 2) - xaux(:, 1);
    f_x = fun(xaux);
    
    z = -delta .* sum((f_x(:, 1:end-1) + f_x(:, 2:end)), 2)/2;


        
    % comp = dilog(1 - x);
end