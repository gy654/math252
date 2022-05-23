A = [2 1  1/2 1/4; 1 4 1 1/2; 1/2 1 4 1; 1/4 1/2 1 2];
[R] = CholeskyFact(A)
function [R] = CholeskyFact(A)
    [n, n2] = size(A);
    R = eye(n);
    RS = eye(n);
    RS(1,1)= A(1,1);
    if RS(1,1)<0
        disp('error')
    else
        R(1,1) = sqrt(RS(1,1));
    end
    for i = 2:n
        for j = 1:i-1
            sum = 0;
            for k = 1:j-1
                sum = sum + R(i,k)*R(j,k);
            end
            R(i,j) = (1/R(j,j))*(A(i,j)- sum);
            R(j,i) = R(i,j);
        end
        sum2 = 0;
        for k = 1:i-1
            sum2 = sum2 +R(i,k)^2;
        end
        RS(i,i)= A(i,i)- sum2;
        if RS(i,i)<0
            disp('error')
        else
            R(i,i)= sqrt(RS(i,i));
        end
    end
end
        

        
