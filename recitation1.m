%% 
% clc: clear the command window
% clear: clear all variables
% close




%% clc; clear; close all

v = [1,2,3];
w = [1;2;3];

% transpose(conjugate) of a vector
z = v';
x = [i,i];
f = @(a) 2*a;
y = f(x);

%% matrix
clc;clear;close all
v = [1:0.5:5]'
% get vector partially
v(1:3)
size(v);

n=9;
e = ones(n,1);
% sp: sparse; diags
A = spdiags([e -2*e e], -1:1,n,n);
% get matrix partial
full(A)
full(A(:,1:2))

w = A*v

% back slash operation is 
vv=A\w;
vvv=inv(A)*w;

B= eye(9);
C = A\B; % A(-1)*B
D = inv(A)*B
% C is the same as D

%% vector, matrix and function
% index starts with 1!!!!

clc;clear; close all

% function handle: treat x as a variable
f = @(x) sin(x);
a = [1:0.1:10]'
b = f(a);
plot(a,b)
title('sin(x)')
xlabel('x')
ylabel('y')
legend('sinx')

% run section

%% Built-in function in Matlab
clc, clear, close all
% QR factorization
n=9;
e=ones(n,1);
A = spdiags([e -2*e e], -1:1,n,n);
% to see a whole matrix instead of a sparse one, use full(A)
[Q, R] = qr(full(A))

%% for loop and while loop
clc, clear, close all

% for loop: "end!"
for n = [1 2 3]
    2*n
    plot()
    hold on
end

for n = [1: 0.5: 3]
    2*n
end


tic % timer
ii=0; % initial value
idx=1;
y =  ones(11,1); % store sin(ii) into y; allocate it before the loop
while ii< 10
    y(idx)=sin(ii)
    ii = ii + 1;
    idx = idx +1;
end
a=toc

%% 
clc; clear; close all
x=1
if x<0
    a=1
elseif 0<x<2
    a=2
else
    a=3
end

%% plot and more interesting features
f = @(x) sin(x);
g = @(x) log(x);

a = [1:0.1:10]';
b = f(a);
c = g(a);

% two functions in the same plot
% first one
plot(a,b,'b--')
hold on 
plot(a,c,'r*')
% second way
plot(a,b,'-o',a,c,'--', [0,1])
% '*', '--' or other things on the plot
% color
title('sin(x) vs. log(x)')
xlabel('x')
ylabel('y')
legend('sinx', 'log(x)')
