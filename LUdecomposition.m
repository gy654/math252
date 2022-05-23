
A=[1,2,3;4,5,6;7,8,9];

[L,U]=LU(A)
function [L, U] = LU(A)
    n= size(A);
    L = eye(n);
    U = eye(n);
    U(1, :)= A(1, :);
    for ii = 2:n
        for jj =1:ii-1
            sum = 0;
            for kk = 1: jj-1
                sum = sum+L(ii,kk)*U(kk,jj);
            end
            L(ii,jj) = (1/U(jj,jj))* (A(ii,jj) - sum);
        end
        for jj = ii: n
            sum2 = 0;
            for kk = 1: ii-1
                sum2 = sum2+L(ii,kk)*U(kk,jj);
            end
            U(ii,jj) = A(ii,jj)-sum2;
        end
    end
end
       
          
    
    