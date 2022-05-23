A = [3 0 0; 0 -5 0; 0 0 4];
a = sqrt(3^2 +5^2 +4^2);
N = 100;
v = zeros(3,N);
v(:, 1) = (1/a)*[1;1;1];
e = zeros(N,0);
for k = 1: N-1
    w = A*v(: ,k);
    v(:, k+1) = (1/norm(w))*w;
    e(k+1) = v(:, k+1)' *A*v(:, k+1);
end
e(N)
v(:, N)



%%
A = [-2 1 4; 1 1 1;4 1 -2];
a = sqrt(1+4+1);
N =100;
v = zeros(3,N);
v(:, 1) =(1/a)*[1;2;-1];
e = zeros(N, 0);
for k = 1: N-1
    w = A*v(:, k);
    v (:,k+1) = (1/norm(w))*w;
    e(k+1) = v(:, k+1)' *A * v(:, k+1);
end
v(:, N)
e(N)

%%
%6(d)
% construct matrix A
v0 = 8 * ones(1, 4);
v1 = 4 * ones(1, 3);
v2 = 2 * ones(1, 2);
v3=  1 * ones(1, 1);

A = diag(v0) +diag(v1, -1) + diag(v1, 1)+ diag(v2, -2)+ diag(v2, 2)+diag(v3, 3) + diag(v3, -3);
full(A)

% use the power method to find the largest eigenvalue and its corresponding
% eigenvector
e_vec = zeros(4, 4);
e_val = zeros(4,4);
N =100;
v = zeros(4,N);
v(:, 1) =(1/2)*[1;1;1;1];
e = zeros(N, 0);
for k = 1: N-1
    w = A*v(:, k);
    v (:,k+1) = (1/norm(w))*w;
    e(k+1) = v(:, k+1)' *A * v(:, k+1);
end
fprintf("lambda 1 = %f \n", e(100))
disp("v1: ")
v(:, 100)
e_vec (:, 1) = v(:, 100);
e_val (1, 1) = e(100);
%%
% set v0 orthorgonal to q1to find the second largest eigenvalue and eigenvector

N =30;
v = zeros(4,N);
ort = null(e_vec(:, 1).');
v(:, 1)= ort(:, 1);
e = zeros(N, 0);
for k = 1: N-1
    w = A*v(:, k);
    v (:,k+1) = (1/norm(w))*w;
    e(k+1) = v(:, k+1)' *A * v(:, k+1);
end
fprintf("lambda 2 = %f \n", e(N))
disp("v2: ")
v(:, N)
e_vec (:, 2) = v(:, N);
e_val (2, 2) = e(N);

%% set v0 orthorgonal to q1 and q2 to find the third largest eigenvalue and eigenvector

N =16;
v = zeros(4,N);
ort = null([e_vec(:, 1).'; e_vec(:, 2).']);
v(:, 1)= ort(:, 1);
e = zeros(N, 0);
for k = 1: N-1
    w = A*v(:, k);
    v (:,k+1) = (1/norm(w))*w;
    e(k+1) = v(:, k+1)' *A * v(:, k+1);
end
fprintf("lambda 3 = %f \n", e(N))
disp("v3: ")
v(:, N)
e_vec (:, 3) = v(:, N);
e_val (3, 3) = e(N);


%% set v0 orthorgonal to q1 and q2 and q3 to find the third largest eigenvalue and eigenvector

N = 7;
v = zeros(4,N);
ort = null([e_vec(:, 1).'; e_vec(:, 2).';e_vec(:,3).']);
v(:, 1)= ort(:, 1);
e = zeros(N, 0);
for k = 1: N-1
    w = A*v(:, k);
    v (:,k+1) = (1/norm(w))*w;
    e(k+1) = v(:, k+1)' *A * v(:, k+1);
end
fprintf("lambda 4 = %f \n", e(N))
disp("v4: ")
v(:, N)
e_vec (:, 4) = v(:, N);
e_val (4, 4) = e(N);
X = e_vec * sqrt(e_val) * e_vec.'
X * X'



