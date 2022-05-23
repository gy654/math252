clc
close all

%%
clc
clear all
close all
x = 10 * rand(100,1)-5;
y = 3 * x + 4 * randn(100,1);
data = [x y];
data = data - mean(data);
% mean(data) calculate the mean of the columns, data-mean(data) makes sure
% the data is at the center


plot(data(:, 1), data(:,2),'o')

% PCA line
[U_MATLAB, S_MATLAB, V_MATLAB] = svds(data, 1); % returns 1 largest singular value
t = -0.25:0.1:25;
PCA = (V_MATLAB*t)'; % just the direction of the largest eigenvector of P^TP/ direction of projection
hold all;
plot(PCA(:, 1), PCA(:,2), 'k-')

% LSq line
% Ax= b, QRx=b, Rx = Q'b
A = [x 0*x+1];
b = y;
[Q, R ]= qr(A);
coefs = R\(Q'*b);
line = A* coefs;  % Ax
plot(x, line, 'r-')

%%
close all
theta = 0:0.01:5*pi;
r = 1.5 * theta;
evar = 1;
x = r.*cos(theta) + evar* randn(size(r)); % parametrization + purturbation
y = r.*sin(theta) + evar* randn(size(r)); % parametrization + purturbation
z = rand(size(r)); % z is randomly generated
% x-y plane, the swiss-roll structure is clearest
data = [x' z' y'];% size : 1571*3
data = data - mean(data);
figure(1)
plot3(data(:,1), data(:,2), data(:,3), 'o')
[U, S, V] = svds(data, 2);
PCA = (U*S);

figure(2)
plot(PCA(:,1), PCA(:, 2), 'o')  % coefficient of  projection on the perspective spanned by the two largest eigenvector of P^TP

figure(3)
r = 2; % only cares about the most significant two singular values
projected_data = U(:, 1:r) * S(1:r, 1:r) * (V(:,1:r))'; 
plot3(data(:,1), data(:,2), data(:,3), 'bo');
hold on
plot3(projected_data(:,1), projected_data(:,2), projected_data(:,3), 'ro')

%% set up least square problem to project onto a plane
A =  data(:, 1:2);
b = data(:,3);
[Q, R]= qr(A);
coefs =  R\(Q'*b);
lsq_data = [data(:,1:2) A * coefs]
figure(3)
plot3(lsq_data(:,1), lsq_data(:,2), lsq_data(:,3), 'go')

%%
%%
%(a) table of data
X = [0 0.5 1.0 1.5 2.0 2.5]';
Y = [0 0.2 0.27 0.30 0.32 0.33]';
coef_a = X.^3;
coef_b = X.^2;
coef_c = X;
coef_d = ones(6,1);

% write coefficients into matrix A
A = zeros(6,4);
A(:,1) = coef_a;
A(:,2) = coef_b;
A(:,3) = coef_c;
A(:,4) = coef_d;


% all six points
[Q, R]= qr(A);
coefs =  R\(Q'*Y);
x = linspace(0,3, 100);
y_lsq = coefs(1)* x.^3 +coefs(2) * x.^2 + coefs(3) * x + coefs(4);
plot(x, y_lsq,'blue')
hold on


% only the data i = 0,1,2,3,4
[Q2, R2]= qr(A(1:5,1:4));
coefs_2 =  R2\(Q2'*Y(1:5));
x = linspace(0,3, 100);
y2_lsq = coefs_2(1)* x.^3 +coefs_2(2) * x.^2 + coefs_2(3) * x + coefs_2(4);
plot(x, y2_lsq)
hold on



% only the data i = 0,1,2,3
[Q3, R3]= qr(A(1:4,1:4));
coefs_3 =  R3\(Q3'*Y(1:4));
x = linspace(0,3, 100);
y3_lsq = coefs_3(1)* x.^3 +coefs_3(2) * x.^2 + coefs_3(3) * x + coefs_3(4);
plot(x, y3_lsq)
hold on
scatter(X, Y, 'filled')
xlabel('x')
ylabel('y')
legend('6 points interpolation','5 points interpolation', '4 points interpolation', 'data points')

% case three is a sqecial least square problem because the least square
% polynomial we find fits perfectly with the actual data point without any error. 

%%
%(b)

clc
clear all
close all
X = [0 0.5 1.0 1.5 2.0 2.5]';
Y = [0 0.2 0.27 0.30 0.32 0.33]';
coef_a = X.^5;
coef_b = X.^4;
coef_c = X.^3;
coef_d = X.^2;
coef_e = X;
coef_f = ones(6,1);

A = zeros(6,6);
A(:,1) = coef_a;
A(:,2) = coef_b;
A(:,3) = coef_c;
A(:,4) = coef_d;
A(:,5) = coef_e;
A(:,6) = coef_f;

[Q, R]= qr(A);
coefs =  R\(Q'*Y);
x = linspace(0,3, 100);
y_lsq = coefs(1)* x.^5 +coefs(2) * x.^4 + coefs(3) * x.^3 + coefs(4)* x.^2 + coefs(5)*x + coefs(6);
plot(x, y_lsq,'blue')
hold on
scatter(X, Y, 'filled')
legend('five degree interpolation', 'data points')


% we would have to use polynomials with degree greater or equal to 5 so
% that the solution interpolates all six data points


%%  4. polynomial interpolation and error estimation
%(a)
clc
clear all
close all
X = [0 1/2 1]';
Y = exp(3*X);
A = zeros(3,3);
A (:, 1)= [1 1 1]';
A (:, 2)= X;
A (:, 3)= X.^2;
coef = linsolve(A, Y);
x = linspace(0,1,100);
P2 = coef(1) + coef(2) * x + coef(3) * x.^2;
f = exp(3*x);
plot(x, P2)
hold on
plot(x, f)

legend('numerical interpolation using monomial basis','f(x)= exp(3x)', 'location', 'best')
title('numerical interpolation of three nodes x0=0, x1=1/2, x2=1 using monomial basis')
xlabel('x');
ylabel('y');




%% (b)

clc
clear all
close all

x = linspace(0,1,100);
L0 = 2* x.^2 - 3*x +1;
L1 = -4*x.^2 + 4*x;
L2 = 2*x.^2-x;

plot(x, L0)
hold on
plot(x, L1)
hold on
plot(x, L2)
title('lagrange basis')
legend('lagrange basis L0', 'lagrange basis L1', 'lagrange basis L2')

%%


clc
clear all
close all

c0 = 1;
c1 = 4 * e^(3/2) - e^3 - 3;
c2 = 2 - 4 * e^(3/2) + 2 * e^3;

x = linspace(0,1,100);
P2_a = c0 + c1 * x + c2 * x.^2;
f = exp(3 * x);

plot(x, f)
hold on
plot(x, P2_a)

title('interpolate f(x) = exp(3x) using three nodes x0=0, x1=1/2, x2=1')
xlabel('x')
ylabel('y')
legend('p2 interpolation', 'f=exp(3x)','location', 'best')


%%
%3(c)
K2 = [];
n1 = 5;
x_5 = linspace(-1,1,n1+1);
A_5 = fliplr(vander(x_5));
K2(1) = cond(A_5);
disp('when n=5: ')
K2(1)

n2 = 10;
x_10 = linspace(-1,1,n2+1);
A_10 = fliplr(vander(x_10));
K2(2) = cond(A_10);
disp('when n=10: ')
K2(2)


n3 = 20;
x_20 = linspace(-1,1,n3+1);
A_20 = fliplr(vander(x_20));
K2(3) = cond(A_20);
disp('when n=20: ')
K2(3)

n4 = 30;
x_30 = linspace(-1,1,n4+1);
A_30 = fliplr(vander(x_30));
K2(4) = cond(A_30);
disp('when n=30: ')
K2(4)

%%
clear; close all; clc
% generate Chebyshev points on the domain [-1,1]
a = -1;
b = 1;
% define the indication function
f = @(x) (x >=0);

max_error = [];
L2_error = [];
n_list = [2 4 8 16 32 64 128 256];
for n = [2 4 8 16 32 64 128 256];
    i = linspace(0,n, n+1);
    x = 0.5 * (a+b) +0.5 * (b-a)*cos((i + 0.5 )* pi/(n+1));
    y = f(x);
    % evaluate error pn-f at a large number of uniformly distributed points (at ~ 10n points)
    x0 = linspace(-1,1,10*n);
    y0 = zeros(length(x0), 1);
    for N = 1:length(x0);
        y0(N) = lagrange_interpolant(x, y, x0(N));
    end
    error = abs(y0 - f(x0)');
    sum_error_sqr = sum(error.^2);
    max_error = [max_error max(error)];
    L2_error = [L2_error sqrt((1/5*n)* sum_error_sqr)];

end
disp('table of maximum error:')
table_max_error = [n_list; max_error]
% expect convergence in the maximum norm

disp('table of L2 error:')
table_l2_error = [n_list; L2_error]
% don't expect convergence in the L2 norm



%%
function y0 = lagrange_interpolant(x, y, x0)
% x i sthe vector of abscissas
% y is the matching vector of ordinates
% x0 represents the target to be interpolated
% y0 represents the solution from the lagrange interpolation
y0 = 0;
n = length(x);
for j = 1:n
    t=1;
    for i = 1: n
        if i ~=j
            t=t * (x0-x(i))/(x(j)-x(i));
        end
    end
    y0 = y0+ t* y(j);
end
end
