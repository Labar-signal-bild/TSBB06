%% Förberedelser

n = 2*rand(3,1)-1;
n = n/norm(n);
alpha = 2*pi*rand(1,1);

Rrod = Rodrigues(n, alpha);

[newn, newalpha] = angleAxis(Rrod);
expR = createExpR(n, alpha); %det(expR) = -1

%rotationmatrix och expR ser ej likadana ut. Det(expR) = -1, Det(rotmat)=1

[quatR, q] = quaternionR(n, alpha);


%%  4.1 

Rrod*n; %Check so that R*n = n

% ANSWER: R*n = n. As expected
%det(Rrod) = 1

n
newn
alpha
newalpha

% NO! n is correct but not alpha HELP

%% 4.2

% ANSWER: No, Rrod and expR is not the same..

[e1 l1]=eig(expR);
[e2 l2]=eig(alpha*liu_crossop(n));

%ANSWER: Yes it's what wwe expected but have to find the right places

%% 4.3

%ANSWER: norm(q) = 1

x0 = randn(3,1);

xprim = sandwichproduct(q, x0);
quatR*x0;

%ANSWER: xprim and quatR*x0 give the same result
%det(quatR) = 1

liu_R_from_q(q);
quatR;

%ANSWER: Yes they are the same


%% 5 Estimation of rigid transformations
 
%% 5.1

N = 9;
x1 = 2*rand(3,N)-1;

t = 2*rand(3,1)-1;      % translation vector
n = 2*rand(3,1)-1;
n = n/norm(n);          % normalized rotation axis
a = 2*pi*rand(1,1);     % rotation angle
R = liu_rodrigues(n,a); % R from (n,a) det(R) = 1

x2 = R*x1+t*ones(1,N);

s = 0.01;
x1n = x1 + s*randn(3,N);
x2n = x2 + s*randn(3,N);

%% 5.2

[Rest, test] = Rsvd(x1n, x2n,0);

[e3 l3] = eig(Rest);

ne3Cross = e3(:,1)'*n; %If zero e3(:,1) is orthogonal with n 

%ANSWEAR: Not a rotation matrix around n!!!!! 

%Rest*n ~= n, hmm är något fel??

x2e = Rest*x1n + test*ones(1,N);
err = norm(x2e - x2n,'fro');

% err = 1.8761. This is not what we expected. The noise we added is of the
% magnitude 0.01 on 9 points. Which would not be summed to 1.8

%% 5.3 A

load rigiddataA

N = 90; %OBS! Max nmb for N == 90!
dataPreA = data1(:,1:N);
dataAftA = data2(:,1:N);

[RdataA, tdataA] = Rsvd(dataPreA,dataAftA,0); 
[ndataA, alphadataA] = angleAxis(RdataA);
% ANSWEAR: det(RdataA) = 1, RdataA*ndataA = ndataA -> SO(3)!

dataAftErrA = RdataA*dataPreA + tdataA*ones(1,N);
dataErrA = norm(dataAftErrA - dataAftA,'fro');

%% 5.3 B

load rigiddataB

N = 90; %OBS! Max nmb for N == 90!
dataPreB = data1(:,1:N);
dataAftB = data2(:,1:N);

[RdataBold, tdataBold] = Rsvd(dataPreB,dataAftB,0); 
[ndataB, alphadataB] = angleAxis(RdataBold);
% ANSWEAR: det(RdataA) = -1, RdataA*ndataA = ndataA -/> SO(3)!

dataAftErrBold = RdataBold*dataPreB + tdataBold*ones(1,N);
dataErrBold = norm(dataAftErrBold - dataAftB,'fro');

%----------------------- AFTER FIXING OOP ---------------------------
[RdataBnew, tdataBnew] = Rsvd(dataPreB,dataAftB,1); 
[ndataB, alphadataB] = angleAxis(RdataBnew);
% ANSWEAR: det(RdataA) = 1, RdataA*ndataA = ndataA -> SO(3)!

dataAftErrBnew = RdataBnew*dataPreB + tdataBnew*ones(1,N);
dataErrBnew = norm(dataAftErrBnew - dataAftB,'fro');
