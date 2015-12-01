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


















