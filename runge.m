% illustration of Runge phenomenen
close all

%f = @(x) 1./(1+(x-2).^2) + 1./(1+(x+2).^2);
%f = @(x) sin(pi*x);
%f = @(x) 1./(1+25*x.^2);
f = @(x) sign(x-0.5);


a = -1;
b = 1;
xVals = a:0.01:b; % Not the point that will be interpolated at
                  % Only for plotting
clf

N = 16; % Maximum polynomial degree that I want to interpolate with
Ns = 10; % Minimum polynomial degree that I want to interpolate with
cols = parula(N-Ns+2);
cc = 0;
leg = cell(1,N-Ns+1);
for n=Ns:N
    cc = cc + 1;
    % interpolation points
    % equidistant
    %x = a+(0:n)*(b-a)/n;
    % Chebyshev
    x = 0.5*(b-a)*cos(pi*(2*(0:n)+1)/(2*n+2));

    coeffs = polyfit(x,f(x),n); %fits polynomial

    fVals = polyval(coeffs,xVals); %evaluates polynomial

    plot(x,f(x),'o','color',cols(cc,:),'markerfacecolor',cols(cc,:),'markersize',10);
    hold all
    pl(cc) = plot(xVals,fVals,'color',cols(cc,:),'LineWidth',2);
    leg{cc} = ['n = ' num2str(n)];
    hold all
    xlabel('x')
    axis([a,b,-2,2])
end
pl(cc+1) = plot(xVals,f(xVals),':k','LineWidth',4);
hold off
leg{cc+1} = 'f(x)';
legend(pl,leg);


%%
syms x;
eqn = 0.694^x == 10^(-2);
S = solve(eqn, x);
