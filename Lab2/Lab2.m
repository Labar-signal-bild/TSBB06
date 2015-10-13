

n = 2*rand(3,1)-1;
n = n/norm(n);
alpha = 2*pi*rand(1,1)

rotationmatrix = Rodrigues(n, alpha);

[newn, newalpha] = angleAxis(rotationmatrix);
newalpha