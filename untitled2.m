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
publish('untitled2.m', 'pdf')