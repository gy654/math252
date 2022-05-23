
U = [1 2 6 -1; 0 3 1 0; 0 0 4 -1; 0 0 0 2];
b = [-1;-3;-2;4];
[x]= backwardSub(U,b)
function [x]  = backwardSub (U, b)
    [n1b,n2b] = size(b);
    [n1U,n2U] = size(U);
    if n2U ~= n1b
        disp('The size of U and b are not compatible')
    else
        n = n1b;
        x = zeros(n, 1);
        x(n) = b(n)/U(n,n);
        for i = 1:n-1
            sum = 0;
            for jj = n-i+1: n
                sum = sum + U(n-i, jj)*x(jj);
            end
            x(n-i)= (b(n-i)-sum)/U(n-i, n-i);
        end
    end
end
   

    
    
    