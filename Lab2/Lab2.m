%% Förberedelser

n = 2*rand(3,1)-1;
n = n/norm(n);
alpha = 2*pi*rand(1,1);

rotationmatrix = Rodrigues(n, alpha);

[newn, newalpha] = angleAxis(rotationmatrix);
expR = createExpR(n, alpha);

%rotationmatrix och expR ser ej likadana ut. Det(expR) = -1, Det(rotmat)=1

[quatR, q] = quaternionR(n, alpha);


%%  4.1 

rotationmatrix*n

% ANSWER: R*n = n. As expected

n;
newn;
alpha;
newalpha;

% NO! n is correct but not alpha

%% 4.2

% ANSWER: No, rotationmatrix and expR is not the same..

[e1 l1]=eig(expR);
[e2 l2]=eig(alpha*liu_crossop(n));

%ANSWER: Yes it's what wwe expected but have to find the right places

%% 4.3

%ANSWER: norm(q) = 1

x0 = randn(3,1);

xprim = sandwichproduct(q, x0);
quatR*x0;

%ANSWER: xprim and quatR*x0 give the same result

liu_R_from_q(q);
quatR;

%ANSWER: Yes they are the same


%% 5 Estimation of rigid transformations
 
%% 5.1

N = 90;
x1 = 2*rand(3,N)-1;

t = 2*rand(3,1)-1;      % translation vector
n = 2*rand(3,1)-1;
n = n/norm(n);          % normalized rotation axis
a = 2*pi*rand(1,1);     % rotation angle
R = liu_rodrigues(n,a) % R from (n,a)

x2 = R*x1+t*ones(1,N);

s = 0.01;
x1n = x1 + s*randn(3,N);
x2n = x2 + s*randn(3,N);

%% 5.2

a0 = mean(x1n,2);
b0 = mean(x2n,2);
A = x1n - a0*ones(1,N); %Räknar ut hut mycket varje punk avviker från mitten av centrioden
B = x2n - b0*ones(1,N); 


[U, S, V] = svd(A*B');
Rest = V*V*U'*U';

%Rest*n ~= n, hmm är något fel??





