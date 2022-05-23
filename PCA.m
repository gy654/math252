x = 10*rand(100,1)-5;
y = 3*x + 4*randn(100,1);
A = [x y];
A = A - mean(A);


[U,S,V] = svd(A);
t = U(:,1)* S(1,1);
proj_point = zeros(2, 100);
vert = zeros(2,100)
for i = 1: 100
    proj_point (:, i)= V(:,1) * t(i);
    vert(:,i)= [x(i);y(i)] - proj_point (:, i);
end


scatter(x,y)
hold on

scatter(proj_point(1,:),proj_point(2, :)) 
%%
theta = 0:0.01:5*pi;
r = 1.5*theta;
evar = 1;
x = r.*cos(theta) + evar*randn(size(r));
y = r.*sin(theta) + evar*randn(size(r));
z = rand(size(r));
data = [x' z' y'];
data = data - mean(data);
plot3(data(:,1),data(:,2),data(:,3),'bo')
A= [x;y;z]';

[U,S,V] = svds(A);
v1=V(:,1);
v2=V(:,2);
score1 = []
for i = 1: size(A)
    score1(i)= A(i,:)'.* v1';
end
