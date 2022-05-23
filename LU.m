%%

% make a grid of points
N = 1e4;
x = linspace(-pi,3*pi,N+1);
x = 1/2*(x(2:end) +x(1:end-1));
dx = x(2)-x(1);


% make the coupling matrix
e = ones(N,1);
D = spdiags([e -2*e e],-1:1,N,N);
D = (1/dx^2)*D;
% periodic boundaries
D(1,end) = (1/dx^2);
D(end,1) = (1/dx^2);

% make the time time stepping matrix
dt = 0.001;
A = speye(N) - 1i*dt*D;

% choose a final time
T_f = 5.0;

%%
%close all
% initial condition
psi = exp(-8*(x-pi).^2)';
tic 
% L, U factorization
[L,U] = lu(A);
for k = 1:floor(T_f/dt)
   %disp(['t = ' num2str(k*dt)])
   % solve the system using LU
   temp = L\psi; %forward sub Ly=b
   psi = U\temp; %backward sub  Ux=y
   % plot the solution
   %plot(x,psi)
   plot3(x,real(psi),imag(psi))
   xlabel('x')
   ylabel('$$Re[\psi]$$')
   zlabel('$$Imag[\psi]$$')
   daspect([1 1 1])
   axis([0 2*pi -1 1 -1 1])
   view([-75 25])
   drawnow
end
toc

%%
psi = exp(-8*(x-pi).^2)';
tic
invA = A\eye(N); 
for k = 1:floor(T_f/dt)
   psi = invA*psi;
   %disp(['t = ' num2str(k*dt)])
end
toc