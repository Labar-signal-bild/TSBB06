%% TSBB06 LabD Filter Optimisation

cd /site/edu/bb/MultidimensionalSignalAnalysis/ComputerExercises/ExerciseD/
addpath /edu/annhj876/Skola/TSBB06/LabD

%% Preparatory

% 1: Even, real = even, real
%    Even, odd = even, complex


%% Optimisation of a 1D filter

x = (-3:3)';
N = 21; % Number of samples. 3 times as many as in x
u = (-(pi - pi/N):(2*pi/N):(pi - pi/N))';

B = exp(-i*u*x'); % Basis matrix. Constructed with the DFT-basis functions in its columns
% Hur väljer man B?
figure(1);
for t=1:size(B,2),
subplot(7,1,t);plot([real(B(:,t)) imag(B(:,t))]);
axis([1 21 -1 1])
end


% Notice that this matrix has one column(7) for each coeficient of the filter, and
% this column is a sampled version of the function exp() for a fixed x

Fideal = (pi/4<abs(u))&(abs(u)<3*pi/4); % The ideal  filter
figure(2);plot(u,Fideal,'o');
title('Fideal')
W = ones(size(u)); % The weighted function (all ones)
figure(3);plot(u,W,'o');
title('Weighted function')


G0 = diag(W)/N;
G = B'*G0*B;

dualc = B'*G0*Fideal; % The dual coordinates of F(u) in the basis B

c = inv(G)*dualc % The dual coordinates are transformed into proper coordinates
figure(4);plot(x,real(c),'o'); % These  are  the  coordinates  of F(u) 
title('Proper coordinates')

F = real(B*c); % F is given as a linear combination between the optimised coefcients and the basis function
figure(5);plot(u,[F Fideal]);
title('F, Fideal')


epsilon = sqrt((Fideal - F)'*G0*(Fideal- F)) % The error of the optimised filter i


% ANSWER: epsilon = 0.2034

%% 3.1 Changing the sampling in the frequency domain

x = (-3:3)';
epsilon = [];
Nrange = 7:2:501;
for N = Nrange;
    u = (-(pi - pi/N):(2*pi/N):(pi - pi/N))';
    B = exp(-i*u*x');
    Fideal = (pi/4<abs(u))&(abs(u)<3*pi/4);
    W = ones(size(u));
    G0 = diag(W)/N;
    G = B'*G0*B;
    dualc = B'*G0*Fideal;
    c1 = inv(G)*dualc;
    F1 = B*c1;
    epsilon = [epsilon sqrt(abs((Fideal - F1)'*G0*(Fideal- F1)))];
end
figure(6);plot(Nrange, epsilon);

% The last figure shows the error as a function of the number of samples in
% the frequency domain, N.

% ANSWER: The vaule that the plot converges to in figure (6) is the
% how much the F and Fideal differs

% ANSWER: When just using 7 samples F and Fideal look exactly the same,
% hence the error is 0. BUT the ideal filter does not look like we want it
% to.

% ANSWER: The reasonable value for the samlpes are 40. Then there is not
% improvement when useing a higer value. 

[real(c) real(c1)]

% ANSWER: Yes there is a big difference in c and c1, at least in varannat
% sample.


%% 3.2 Changing the sampling in the signal domain

epsilon = []; %Stores the error for each iteration
M = [];    %Stores the number of coefficients for each iteration
cs = cell(0); %Stores the optimized coefficients for each iteration
figure(7);
for P = 3:11;
    M = [M 2*P+1];  %M is the number of filter coefficients
    N = 10*(2*P+1)+1;     %N is the number of frequency samples
    u = (-(pi - pi/N):(2*pi/N):(pi - pi/N))';
    Fideal = (pi/4<abs(u))&(abs(u)<3*pi/4);
    W = ones(size(u));
    G0 = diag(W)/N;
    x = (-P:P)';
    B = exp(-i*u*x');
    G = B'*G0*B;
    dualc = B'*G0*Fideal;
    c = inv(G)*dualc;
    cs{end+1} = c;
    F = B*c;
    epsilon = [epsilon sqrt(abs((Fideal - F)'*G0*(Fideal- F)))];
    subplot(3,3,P-2);plot(u,real(F));
end
figure(8);plot(M,epsilon);

% Figure (8) shows how the error depends on the number of filter coeficients.

% ANSWER: When increasing x we get more basises hence we get more övertoner 
% and that gives us a better approximation

% ANSWER: All the övertoner does not exist in the filter we're trying to
% approximate

%% 3.3 Changing the frequency weighting function

x = (-11:11)';
N = 231; % Number of samples. 3 times as many as in x
u = (-(pi - pi/N):(2*pi/N):(pi - pi/N))';

B = exp(-i*u*x'); % Basis matrix. Constructed with the DFT-basis functions in its columns

Fideal = (pi/4<abs(u))&(abs(u)<3*pi/4); % The ideal  filter
figure(2);plot(u,Fideal,'o');
title('Fideal')

a1 = exp(-(u+1.5).^10);
a2 = a1(end:-1:1);
W1 = (a1+a2);

W2 = 1-(a1+a2);

a3 = exp(-3*abs(u));
W3 = a3;


[F1 c1 G01] = getF(W1,Fideal,N,B);
[F2 c2 G02] = getF(W2,Fideal,N,B);
[F3 c G0] = getF(W3,Fideal,N,B);

figure(11);
subplot(2,1,1);
plot(u,W1,'o');
title('Weighted function 1')
subplot(2,1,2);
plot(u,[F1 Fideal]);

figure(12);
subplot(2,1,1);
plot(u,W2,'o');
title('Weighted function 2')
subplot(2,1,2);
plot(u,[F2 Fideal])

figure(13);
subplot(2,1,1);
plot(u,W3,'o');
title('Weighted function 3')
subplot(2,1,2);
plot(u,[F3 Fideal])

epsilon = sqrt(abs((Fideal - F3)'*G0*(Fideal- F3)))
% epsilon = 0.0288

% F does not look like we expcted in all cases


%% Allmän skit

x=[-3:0.1:3];
a = exp(-(u+1.5).^10);
a2 = a(end:-1:1);
gauss_a = 1-(a+a2);
figure(10)
plot(u, gauss_a)

% How does this influence F? Don't thing it gets better/worse...

a = exp(-3*abs(u));

%% 3.4 Spatial mask

figure(22);plot(x,real(c),'o');grid

mask = zeros(1,23)';
A = [2, 6, 10, 12, 14, 18, 22];
mask(A) = 1;
S = mask2matrix(mask);

% How is the matrix S constructed from the mask? 
% ANSWER: One 1 on each row, the rest on that row is 0.


Bmask = B*S;


Gmask = Bmask'*G0*Bmask;
dualc = Bmask'*G0*Fideal;
c = inv(Gmask)*dualc;
F = Bmask*c;
epsilon = sqrt(abs((Fideal - F)'*G0*(Fideal- F)))

% epsilon = 0.0296
% epsilondiff = 8.0000e-04. It has increased since we took some of the
% weights away

figure(24);plot(u,[F Fideal]);
title('F, Fideal')

figure(23);plot(x,S*real(c),'o');grid

%% The same code as in 3.1 but different W

x = (-3:3)';
N = 21; % Number of samples. 3 times as many as in x
u = (-(pi - pi/N):(2*pi/N):(pi - pi/N))';

B = exp(-i*u*x'); % Basis matrix. Constructed with the DFT-basis functions in its columns
% Hur väljer man B?

Fideal = (pi/4<abs(u))&(abs(u)<3*pi/4); % The ideal  filter
figure(2);plot(u,Fideal,'o');
title('Fideal')
W = W3(116-10:116+10); % The weighted function (all ones)
figure(3);plot(u,W,'o');
title('Weighted function')


G0 = diag(W)/N;
G = B'*G0*B;

dualc = B'*G0*Fideal; % The dual coordinates of F(u) in the basis B

c = inv(G)*dualc % The dual coordinates are transformed into proper coordinates
figure(4);plot(x,real(c),'o'); % These  are  the  coordinates  of F(u) 
title('Proper coordinates')

F = real(B*c); % F is given as a linear combination between the optimised coefcients and the basis function
figure(5);plot(u,[F Fideal]);
title('F, Fideal')


epsilon = sqrt((Fideal - F)'*G0*(Fideal- F)) % The error of the optimised filter i


% ANSWER: epsilon =  0.1665

% This is worse then with the spacial mask because at the interval -3:3 W3
% wants to do all the samples good and nothing gets good.


%% 4 Optimisation of a 2D filter

x = (-3:3);
x1 = ones(7,1)*x;
x2 = x1';
N = 51;
u = -(pi-pi/N):(2*pi/N):(pi-pi/N);
u1 = ones(N,1)*u;
u2 = u1';
x1 = x1(:);x2 = x2(:); %Here matrices are reshaped into columns
u1 = u1(:);u2 = u2(:);

s = 1;
Fideal=exp(-((u1-pi/2).^2+(u2-pi/2).^2)/s^2)+...
exp(-((u1+pi/2).^2+(u2+pi/2).^2)/s^2);
r = sqrt(u1.^2+u2.^2);
W = (r + 0.1).^(-1);
G0 = diag(W);
B = exp(-i*(u1*x1'+u2*x2'));
G = B'*G0*B;
dualc = B'*G0*Fideal;
c = inv(G)'*dualc;
F = real(B*c);     %Since imag = 0
epsilon = sqrt(abs((Fideal-F)'*G0*(Fideal-F)))

figure(14);
subplot(3,1,1);mesh(u,u,reshape(Fideal,N,N));
subplot(3,1,2);mesh(u,u,reshape(F,N,N));
subplot(3,1,3);mesh(u,u,reshape(abs(Fideal-F),N,N));
c2D=reshape(real(c),7,7)
figure(15);mesh(x,x,c2D);


%%
clear all

x = (-2:2);
x1 = ones(5,1)*x;
x2 = x1';
N = 51;
u = -(pi-pi/N):(2*pi/N):(pi-pi/N);
u1 = ones(N,1)*u;
u2 = u1';
x1 = x1(:);x2 = x2(:); %Here matrices are reshaped into columns
u1 = u1(:);u2 = u2(:);

s = 1;
Fideal=exp(-((u1-pi/2).^2+(u2-pi/2).^2)/s^2)+...
exp(-((u1+pi/2).^2+(u2+pi/2).^2)/s^2);
r = sqrt(u1.^2+u2.^2);
W = (r + 0.1).^(-1);
G0 = diag(W);
B = exp(-i*(u1*x1'+u2*x2'));
G = B'*G0*B;
dualc = B'*G0*Fideal;
c = inv(G)'*dualc;
F = real(B*c);     %Since imag = 0
epsilon = sqrt(abs((Fideal-F)'*G0*(Fideal-F)))

figure(16);
subplot(3,1,1);mesh(u,u,reshape(Fideal,N,N));
subplot(3,1,2);mesh(u,u,reshape(F,N,N));
subplot(3,1,3);mesh(u,u,reshape(abs(Fideal-F),N,N));
c2D=reshape(real(c),5,5)
figure(17);mesh(x,x,c2D);





