function muest = newton_rapson(lambda_trad, c)
    M = size(lambda_trad, 1);
    N = M / c;
    p = 1;
    for ii = 1:M
        p = conv(p,[-1,lambda_trad(ii)]);
    end
    for m = 1:M
        ptemp = 1;
        for l = 1:M
            if l~=m
                ptemp = conv(ptemp,[-1,lambda_trad(l)]);
            end
        end
        p = p-[0,1/N*lambda_trad(m)*ptemp];
    end
    p = [p,0]; 
    mu_init_l = sort(roots(p),'ascend');


    err=1e-6;
    M = size(lambda_trad, 1);
    N = M / c;

    muest=zeros(M, 1);

    for m=1:M,
        if m>1,
            mutemp=(lambda_trad(m)+lambda_trad(m-1))/2;
            mutemp2=(2*lambda_trad(m)+lambda_trad(m-1))/3;
            while (abs((mutemp-mutemp2)./mutemp)>err) || (abs(sum(lambda_trad./(lambda_trad - mutemp)) - N ) > 1e-3),
                mutemp2=mutemp;
                ftemp=sum(lambda_trad./(lambda_trad-mutemp))/M-1/c;
                fdtemp=sum(lambda_trad./((lambda_trad-mutemp).^2))/M;
                mutemp=mutemp-ftemp/fdtemp;
                ind=1;
                while mutemp>lambda_trad(m) | mutemp<lambda_trad(m-1),
                    mutemp=mutemp2-1/ind*ftemp/fdtemp;
                    ind=ind+1;
                end;
            end;
        else,
            % mutemp=lambda_trad(m)/2;
            mutemp=mu_init_l(m);
            mutemp2=2*lambda_trad(m)/3;
            while abs((mutemp-mutemp2)./mutemp)>err || (abs(sum(lambda_trad./(lambda_trad - mutemp)) - N ) > 1e-3),
                mutemp2=mutemp;
                ftemp=sum(lambda_trad./(lambda_trad-mutemp))/M-1/c;
                fdtemp=sum(lambda_trad./((lambda_trad-mutemp).^2))/M;
                mutemp=mutemp - ftemp/fdtemp;
                ind=1;
                while mutemp>lambda_trad(m) && ind < 1e5,
                    mutemp=mutemp2 - (1/ind)*ftemp/fdtemp;
                    ind = ind + 1;
                end;
            end;
        end;
        muest(m)=mutemp;
    end;
end