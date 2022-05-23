%% Week 1 Notes
% download Matlab newest version, free from NYU account 
% ".m" --> new script # mainly used


% Good Habits:
% Add this to at the beginning of your projects
clc; clear; close all
% clc: clear the command window
% clear: claer all the variables
% close all: close the figures

% create a new section and run it: 
% use "%%", run section if you only want to debug one section 
% use "help" in the command window


% Next Week:
% how to write a function 
% how to publish ".m" file
% how to debug your projects

% Don't get too frustrated for the first time! 
% Things will be better if you are familar with it!

%% Vector, matrix and f(x)

% !!!!!!!!!!! vector starts with index 1, not 0 !!!!!!!!!!!!!

clc; clear; close all
% ";" at the end: the result will not show up in the command window

% row vector
v = [1,2,3];
% col vector
w = [1;2;3];

% conjugate transpose of a vector: v'
z = v';

% i is the imaginary number
x = [i, -i];
y = x';

f = @(a) 2*a;
y = f(x);


%% matrix 
clc; clear; close all

% 1 to 5 with interval 0.5
v = [1:0.5:5]'; 
% get the vector partially 
v(1:3);
size(v);

n = 9;
e = ones(n,1);

% sp: sparse matrix; diags: 
A = spdiags([e -2*e e],-1:1,n,n);
full(A);
w = A*v;

% get the matrix partially 
full(A(:, 1:2));

% "\" is matrix\vector or matrix 
vv = A\w;          % A^(-1)*w
vvv = inv(A)*w;

B = eye(9);
C = A\B;   % A^(-1)*B
D = inv(A)*B;

%% f(x) 
clc; clear; close all

% function handle: treat x as a variable
f = @(x) sin(x);

a = [1:0.1:10]';
b = f(a);
plot(a, b)

title('sin(x)')
xlabel('x')
ylabel('y')
legend('sinx')


%% Built-in funciton in Matlab
clc; clear; close all

% QR Factorization 
n = 9;
e = ones(n,1);
A = spdiags([e -2*e e],-1:1,n,n);

% To see the whole matrix instead of a sparse one, use "full(A)"
[Q,R] = qr(full(A));

%% For Loop and While Loop 
clc; clear; close all

% for loop: "end"!
for n = [1:0.5:3]
    2*n;
end

tic % timer in Matlab
ii = 0; % initial value 
idx = 1; 
y = ones(11,1); % store sin(ii) in to y; allocate it before the loop 
while ii <= 1000
    y(idx) = sin(ii); 
    ii = ii + 1;
    idx = idx + 1;
end
t = toc % end of the timer

%% IF ELSE Statement
clc; clear; close all
x = 1;
if x < 0 
    a = 1
elseif 0 < x < 2
    a = 2
else
    a = 3
end

%% Plot and More interesting features
clc; clear; close all

f = @(x) sin(x);
g = @(x) sqrt(x); 

a = [1:0.1:10]';
b = f(a);
c = g(a);

% plot two function on the same plot:
% 1. plot(x1,y1,x2,y2)
% 2. hold on, hold off --- useful if you plot with For Loop!
plot(a, b, 'o', a, c, '--')
% hold on 
% plot(a, c, 'r*')

% plot on an interval: help plot 

% Plot Features:
% '*', '--' or other things on the plot 
% color 

title('sin(x) vs. sqrt(x)')
xlabel('x')
ylabel('y')
legend('sinx', 'sqrt(x)')










