clc; clear; close all

%% Problem1 (a)
close all;
clear;
clc;
f=@(x) x^3-3*x^2+3; %the function
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
%% Problem1 (b)
close all;
clear;
clc;
x(1) = 1.5;
ff = @(xx) xx^3-3*xx^2+3;
ffprime = @(xx) 3*xx^2 - 6*xx;
for i = 1:6
    % The actual Newton step
    x(i+1) = x(i) - ff(x(i))/ffprime(x(i));
    fprintf('x(%d)=%.16f\n', i-1,x(i));
end
%% Problem1 (c)
close all;
clear;
clc;
x(1) = 2.1;
ff = @(xx) xx^3 -3*xx^2+3;
ffprime = @(xx) 3*xx^2 - 6*xx;
for i = 1:8
    % The actual Newton step
    x(i+1) = x(i) - ff(x(i))/ffprime(x(i));
    fprintf('x(%d)=%.16f\n', i-1,x(i));
end


%% problem 4
close all;
clear;
clc;
x(1) = 3;
ff = @(xx) (xx/3)^2+(1-0.15*xx)^2-1;
ffprime = @(xx) (2/3)*(xx/3)+(2-0.3*xx)*(-0.3/2);
tol = 1e-6;
disp("x( 1) =3")
for i = 1:100
    % The actual Newton step
    x(i+1) = x(i) - ff(x(i))/ffprime(x(i));
    fprintf('x(%2d) =%.16f \n',i+1, x(i+1));
    if (abs(ff(x(i+1)))<tol)
        fprintf('Converged after %d iterations.\n', i);
        fprintf('x(%2d) =%.16f \n',i+1, x(i+1));
        fprintf('Rt = [%.16f, %.16f]',x(i+1), 2-0.3*x(i+1));
        break
    end
    
end

%% problem 5
close all;
clear;
clc;
x(1) = 0;
y(1) = 2;
u(1) = 1;
v(1) = -0.3; % set starting point and initial ray direction
tol = 1e-6; 
for ii = 1:50
    ff = @(tt) ((x(ii)+tt*u(ii))/3)^2+((y(ii)+tt*v(ii))/2)^2-1;
    ffprime = @(tt) 2*u(ii)*(x(ii)+tt*u(ii))/9 + (y(ii)+tt*v(ii))*v(ii)/2;
    t(1)=2;
    for jj = 1:100   
        t(jj+1)= t(jj) - ff(t(jj))/ffprime(t(jj));
        if (abs(ff(t(jj+1)))<tol)
            % fprintf('t(%3d) =%.16f \n',jj+1, t(jj+1));
            % fprintf('Rt = [%.16f, %.16f]\n',x(ii)+u(ii)*t(jj+1), y(ii)+v(ii)*t(jj+1));
            break
        end
    end   
    x(ii+1) = x(ii)+t(jj+1)*u(ii);
    y(ii+1) = y(ii)+t(jj+1)*v(ii);
    n1(ii+1) = (2/3)*x(ii+1);% first component of the normal vector
    n2(ii+1) = (3/2)*y(ii+1);% second component of the normal vector
    N = [n1(ii+1), n2(ii+1)];
    unit_n1 = (2/3)*x(ii+1)/(((2/3)*x(ii+1))^2+((3/2)*y(ii+1))^2)^(1/2);  % first component of the unit normal vector   
    unit_n2 = (3/2)*y(ii+1)/(((2/3)*x(ii+1))^2+((3/2)*y(ii+1))^2)^(1/2);  % second component of the unit normal vector
    unitN = [unit_n1, unit_n2];
    w1(ii) = u(ii)/sqrt((u(ii))^2+(v(ii))^2);
    w2(ii) = v(ii)/sqrt((u(ii))^2+(v(ii))^2);
    w = [w1(ii), w2(ii)];
    V = w- 2 * dot(w,unitN)*(unitN); % ray direction of the next iteration
    u(ii+1) = V(1);
    v(ii+1) = V(2);
end
P = [x', y'] % 50 pairs of points on the ellipse

% plot the ellipse
a=3;
b=2; 
x0=0; 
y0=0;
t=-pi:0.01:pi;
x=x0+3*cos(t);
y=y0+2*sin(t);
plot(x,y,'lineWidth', 3 );
hold on;

% plot 49 ray tracing lines
for ii = 1:49
    plot([P(ii,1), P(ii+1,1)], [P(ii, 2), P(ii+1, 2)], 'r');
    hold on;
end
title('Ray tracing algorithm with 50 steps');
xlabel('x');
ylabel('y');
hold off;

%% problem 6 (a)
close all;
clear;
clc;
ff = @(x) x^5-3*x^2+1;
A = -2:0.25:2; % the interval where we want to find the root. note that roots are at least 0.25 away from each other
tol = 1e-6;
for ii = 1:16 % 17 subintervals for bisection method altogether
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
%% Problem6 (b)
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


