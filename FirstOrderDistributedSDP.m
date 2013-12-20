clear all
close all
clc

Nagent  = 2;
Ncouple = 1;

N = 3;

for k = 1:Nagent
    A(:,:,k) = random('norm',0,1,N,N);
    a(k)     = random('norm',0,1);
end

for k = 1:Nagent
    C(:,:,k) = random('norm',0,1,N,N);
end

lambda = zeros(Ncouple,1);

for iter = 1:50
    for i = 1:Nagent
        cvx_begin

            variable X(N,N);

            minimize( norm_nuc(X) + lambda*trace(C(:,:,i)*X) );
            subject to
                trace(A(:,:,i)*X) == a(i);
                X == semidefinite(N);

        cvx_end
        Xall(:,:,i) = X;
    end
    res = trace(C(:,:,1)*Xall(:,:,1)+C(:,:,2)*Xall(:,:,2));
    lambda = lambda + 1e-1*res;
    res_norm(iter) = norm(res);
end
semilogy(res_norm)