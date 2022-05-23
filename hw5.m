%% 1
M = [pi pi^2/2 pi^3/3;
    pi^2/2 pi^3/3 pi^4/4;
    pi^3/3 pi^4/4 pi^5/5;];
b = [2 pi pi^2-4]';
c = linsolve(M, b);
x = linspace(0, pi, 100);
P2_l2 = c(1) + c(2) * x + c(3) * x.^2;
P2_int = (-4/pi^2) * (x.^2 - pi *x);
f = sin(x);
plot(x, P2_l2)
hold on
plot(x, P2_int)
hold on
plot(x, f)
legend('quadratic L2 inner product','quadratic interpolation','f(x)=sin(x)')
title('Quadratic approximation of sine')
xlabel('x')
ylabel('y')


%% 4(a)
% write a code to compute the definite integral of a function using the
% composite Simpson rule with a given number of points n
clc;
clear all;
close all;
f = @(x) exp(x).* exp(x) .* sin(x) .* sin(x);
a = 0;
b = 1;
n = 6; % number of intervals (must be even)
y = simpson(f,n,a,b)

% error in Simpson's rule: E2(f) = (1/90) * (b-a/2)^5 * (M4/n^4)
% error <  1.2688 - 1.2674 = 0.0014
% when n = 2
error = abs(1.2668 - simpson(f,2,a,b)) % 0.0217
% the error decays at the rate of n^4
% i.e when n is multiplied by k, the error is multiplied by (1/k^4)
desired_error = abs(1.2668 - 1.2674);
k = (error/desired_error)^(1/4);
num_int = 2 * ceil(k);
sprintf('%d points are needed to get 4 decimal places correctly',num_int+1)
% number of points = number of interval +1


%% 4(b)

w = [5/9 8/9 5/9];
x = [-sqrt(3/5) 0 sqrt(3/5)];
f = @(x) exp((x+1)/2).* exp((x+1)/2) .* sin((x+1)/2) .* sin((x+1)/2);
I = (1/2)*(w(1)* f(x(1)) + w(2)* f(x(2)) + w(3)* f(x(3)))

% The result of 3-point Gaussian quadrature rule is 1.2670, so we get 3 correct decimal digits.




%% 5(b)
clc;
clear all;
close all;
p = [1 -1.5 0.6 -1/20];
x = roots(p)
fun1 = @(x) (x-0.5).*(x-0.8873)/((0.1127-0.5)*(0.1127 - 0.8873));
fun2 = @(x) (x-0.1127).*(x-0.8873)/((0.5-0.1127)*(0.5 - 0.8873));
fun3 = @(x) (x-0.1127).*(x-0.5)/((0.8873-0.1127)*(0.8873 - 0.5));

xmin = 0;
xmax = 1;
w0 = integral(fun1,xmin,xmax)
w1 = integral(fun2,xmin,xmax)
w2 = integral(fun3,xmin,xmax)


%% 5(c)
% simpson's rule
y_simp = [];
a = 0;
b = 1;
n = 2;
exact = [1 5/6 3/4 7/10 2/3 9/14 5/8];
k = [1 2 3 4 5 6 7]
for i = [1 2 3 4 5 6 7]   
    f = @(x) x.^i + x;
    I_simp(i) = simpson(f, n, a, b);
end
error_simp = abs(exact - I_simp)
loglog(k, error_simp)
xlabel('x')
ylabel('y')
title("Simpson's rule")

% Gaussian quadrature rule

x = [0.1127 0.5 0.8873];
w = [0.2778 0.4444 0.2778];
exact = [1 5/6 3/4 7/10 2/3 9/14 5/8];
k = [1 2 3 4 5 6 7];
I_gau = [];

for i = [1 2 3 4 5 6 7] 
    f_x =  x.^i + x;  
    I_gau(i) = f_x(1) * w(1) + f_x(2) * w(2) + f_x(3) * w(3);
end
error_gau = abs(exact - I_gau)

%loglog(k, error_gau)
%xlabel('x')
%ylabel('y')
%title('Gaussian quadrature rule')

%% 6(a)
clc;
clear all;
close all;
a = 0.1;
b = 1;
f = @(x) sqrt(x);
m = [10 20 40 80 160];
exact = (2/3) - (1/(15*sqrt(10)));
error_tra = [];
error_sim = [];
for i = [10 20 40 80 160]
    trap_I = trapz( f, a, b, i);
    simp_I = simpson(f, i, a, b);
    error_tra = [error_tra abs(trap_I - exact)];
    error_sim = [error_sim abs(simp_I - exact)];
end
loglog(m, error_tra)
xlabel('m')
ylabel('error')
title('trapezoid rule error')

%%
loglog(m, error_sim)
xlabel('m')
ylabel('error')
title('simpson rule error')


% 6(b)
% composite trapezoidal rule
log_error_tra = log(error_tra);
log_m = log(m);
A = zeros(length(m), 2);
A(:,1) = log_m;
A(:,2) = ones([1, length(m)]);
[Q, R ]= qr(A);
coefs_tra = R\(Q'* log_error_tra');
k_tra = coefs_tra(1)
D_tra = coefs_tra(2)


% composite simpson rule
log_error_sim = log(error_sim);
coefs_sim = R\(Q'* log_error_sim');
k_sim = coefs_sim(1)
D_sim = coefs_sim(2)

%% 6(c) 
% repeat using a =0 instead of a = 0.1
clc;
clear all;
close all;
a = 0;
b = 1;
f = @(x) sqrt(x);
m = [10 20 40 80 160];
exact = integral(f, a, b);
error_tra = [];
error_sim = [];
for i = [10 20 40 80 160]
    trap_I = trapz( f, a, b, i);
    simp_I = simpson(f, i, a, b);
    error_tra = [error_tra abs(trap_I - exact)];
    error_sim = [error_sim abs(simp_I - exact)];
end

loglog(m, error_tra)
xlabel('m')
ylabel('error')
title('trapezoid rule error')

%%
loglog(m, error_sim)
xlabel('m')
ylabel('error')
title('simpson rule error')

% composite trapezoidal rule
log_error_tra = log(error_tra);
log_m = log(m);
A = zeros(length(m), 2);
A(:,1) = log_m;
A(:,2) = ones([1, length(m)]);
[Q, R ]= qr(A);
coefs_tra = R\(Q'* log_error_tra');
k_tra = coefs_tra(1)
D_tra = coefs_tra(2)


% composite simpson rule
log_error_sim = log(error_sim);
coefs_sim = R\(Q'* log_error_sim');
k_sim = coefs_sim(1)
D_sim = coefs_sim(2)


function y = trapz(f, a, b, m)
x = linspace(a, b , m+1);
h = (b-a)/m;
x1 = x(2:m);
sum1 = sum(f(x1));
y = (h/2) * (f(x(1))+f(x(m+1))) + h * sum1;
end

function y = simpson(f, n, a, b)
x = linspace(a, b, n+1);
h = (b-a)/n;
x4 = linspace(x(2), x(n), n/2);
x2 = linspace(x(3), x(n-1), (n-2)/2);
sum4 = 4 * sum(f(x4));
sum2 = 2 * sum(f(x2));
y = (h/3) * (f(x(1)) + sum4 + sum2 + f(x(n+1)));
end 


%%
clc; clear; close all;
x = linspace(-10,10,10000);
y = x.^2;
plot(x, y)
