
close all;
clear all;
clc;
f=@(x) x^3-3*x^2+3; %Enter the Function here
x0 = 1;
x1 = 2;
fprintf("x(0)=1\n")
fprintf("x(1)=2\n")
for i=1:4
    f0=f(x0); %Calculating the value of function at x0
    f1=f(x1); %Calculating the value of function at x1
    y=x1-((x1-x0)/(f1-f0))*f1; %[x0,x1] is the interval of the root
    err=abs(y-x1);
    fprintf("x(%d)=%.16f\n ",i+1,y)
    x0=x1;
    x1=y;
end

close all;
clear all;
clc;
x(1) = 1.5;
ff = @(xx) xx^3-3*xx^2+3;
ffprime = @(xx) 3*xx^2 - 6*xx;
for i = 1:6
    % The actual Newton step
    x(i+1) = x(i) - ff(x(i))/ffprime(x(i));
    fprintf('x(%d)=%.16f\n', i-1,x(i));
end

x(1) = 2.1;
ff = @(xx) xx^3 -3*xx^2+3;
ffprime = @(xx) 3*xx^2 - 6*xx;
for i = 1:8
    % The actual Newton step
    x(i+1) = x(i) - ff(x(i))/ffprime(x(i));
    fprintf('x(%d)=%.16f\n', i-1,x(i));
end



x(1) = 3;
ff = @(xx) (xx/3)^2+(1-0.15*xx)^2-1;
ffprime = @(xx) (2/3)*(xx/3)+(2-0.3*xx)*(-0.3/2);
tol = 1e-6;
disp("x( 1) =3")
for i = 1:100
    % The actual Newton step
    x(i+1) = x(i) - ff(x(i))/ffprime(x(i));
    fprintf('x(%2d) =%.16f \n',i+1, x(i+1))
    if (abs(ff(x(i+1)))<tol)
        fprintf('Converged after %d iterations.\n', i)
        fprintf('x(%2d) =%.16f \n',i+1, x(i+1))
        fprintf('Rt = [%.16f, %.16f]',x(i+1), 2-0.3*x(i+1))
        break
    end
    
end



ff = @(x) x^5-3*x^2+1;
A = -2:0.25:2;
tol = 1e-6;
for ii = 1:16
    a = A(ii);
    b = A(ii+1);
    if ff(a)*ff(b)>0
        root(ii) = NaN;
    else
        nits = 0; 
        a_k = a; 
        b_k = b; 
        m_k = (b_k+a_k)/2; 
        while abs(ff(m_k)) > tol 
            nits = nits + 1; 
            m_k = (b_k+a_k)/2; 
            if sign(ff(m_k)) == sign(ff(a_k)) 
                a_k = m_k;
            else
                b_k = m_k;
            end
            break
        end
        root(ii) = m_k; 
    end
end
root=(root(~isnan(root)))
%(b)
f = @(x) x^5-3*x^2+1;
fprime = @(x) 5*x^4-6*x;
x(1) = -0.62;
tol = 1e-12;
for i = 1:100
    x(i+1) = x(i) - f(x(i))/fprime(x(i));
    if abs(x(i+1)-x(i))<tol
        fprintf('Converged after %d iterations.\n', i)
        fprintf('The first root is: %.12f\n', x(i+1))
        break;
    end
end

x(1) = 0.62;
tol = 1e-12;
for i = 1:100
    x(i+1) = x(i) - f(x(i))/fprime(x(i));
    if abs(x(i+1)-x(i))<tol
        fprintf('Converged after %d iterations.\n', i)
        fprintf('The second root is: %.12f\n', x(i+1))
        break;
    end
end

x(1) = 1.37;
tol = 1e-12;
for i = 1:100
    x(i+1) = x(i) - f(x(i))/fprime(x(i));
    if abs(x(i+1)-x(i))<tol
        fprintf('Converged after %d iterations.\n', i)
        fprintf('The third root is: %.12f\n', x(i+1))
        break;
    end
end

publish('assignment1_grace.m', 'pdf')