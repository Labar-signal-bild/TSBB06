%% Lab1 uppgift 4
% SVAR 4: ca 2px nogranhet

%%

img1=imread('img1.ppm');
img2=imread('img2.ppm');

figure(1);imagesc(img1);
figure(2);imagesc(img2);

y1 = vgg_get_homg([33 279; 267 447; 406 207; 528 223; 638 296; 635 6; 629 416; 742 178]');
y2 = vgg_get_homg([77 410; 320 502; 356 254; 450 241; 548 280; 464 39; 576 384; 585 162]');

figure(1);hold('on');plot(y1(1,:),y1(2,:),'go')
figure(2);hold('on');plot(y2(1,:),y2(2,:),'go')

%% 4.1

N = length(y1);

A = [];

for cnt=1:N
    u = [y1(1,cnt) y2(1,cnt)];
    v = [y1(2,cnt) y2(2,cnt)];
    
    A = [A; [u(1) 0 -u(1)*u(2) v(1) 0 -v(1)*u(2) 1 0 -u(2)];
            [0 u(1) -u(1)*v(2) 0 v(1) -v(1)*v(2) 0 1 -v(2)]];
end


%% 4.2 Imhomogeneous solution

%A = A(1:8,:); % Uncomment this to create a H with minimum points.

A0 = A(:,1:8);
b = A(:,9);

z = -(A0'*A0)^-1*A0'*b;
z = [z ;1];


H1min = [z(1:3) z(4:6) z(7:9)];



%%

H1 = H1min;
y2b = vgg_get_nonhomg(H1*y1);
y1b = vgg_get_nonhomg(inv(H1)*y2);
figure(1);plot(y1b(1,:),y1b(2,:),'rx');
figure(2);plot(y2b(1,:),y2b(2,:),'rx');

%ANSWER: The H1min transform for the minimum number of points is correct, by 
%definition the H matrix is correct if A*z=0, in our case: A*z ~= 0 (s175)
%ANSWER: Very close, about half a pixel wrong (looked at the images)
%ANSWER: A*z ~= 0.1 - 0.8

%%  Symmetric geometric error

e1 = 0;
for k = 1:length(y1),
    e1 = e1 + norm(vgg_get_nonhomg(y2(:,k))-y2b(:,k))^2 + ...
        + norm(vgg_get_nonhomg(y1(:,k))-y1b(:,k))^2;
end
e1 = sqrt(e1);

%ANSWER: Expected error: sqrt(sum(abs(A*z))) = 1.9698. e1 = 1.9029 Kind of the same as expected


%% Un-symmetric geometric error


e2 = 0;
for k = 1:length(y1),
    e2 = e2 + norm(vgg_get_nonhomg(y1(:,k))-y1b(:,k))^2;
end

e2 = sqrt(e2);


e3 = 0;
for k = 1:length(y1),
    e3 = e3 + norm(vgg_get_nonhomg(y2(:,k))-y2b(:,k))^2;
end

e3 = sqrt(e3);


%ANSWER: e2 = 1.4213, e3 = 1.2653.
% The error in image1 and image2 are not the same! What was expected from
% prep question 8

%% Transform image 1 in accordance to the estimated homography

img2t=image_resample(double(img1),H1,640,800); %640, 800 size of image1
figure(3);imagesc(uint8(img2t))

%ANSWER: Yes the resaulting image is eaqual to the overlapping region of
%image 1 and 2

%% Difference between image 1 and 2

figure(4);imagesc(uint8(img2t)-img2)

%ANSWER: The lines show where the two images are different. The lines are
%located where theres a colorchange. In the bottom right corner of 
%figure 4 there's a lot of color, that's because in image1 this isn't even shown.

%There is a car in the picture who is moving. So in the second image the
%car is covering more of the picture then in image1

%% 4.3 Homogeneous solution
A = A(1:8,:); %A0

[U S V] = svd(A);
H3 = reshape(V(:,end),3,3);

z3 = [H3(1:3),H3(4:6),H3(7:9)]';

minError = sum(abs(A*z3)); %Corresponding calc for the minimal case
%minError = 0



%% Geometric error as in 4.2

y2b3 = vgg_get_nonhomg(H3*y1);
y1b3 = vgg_get_nonhomg(inv(H3)*y2);

e3 = 0;
for k = 1:length(y1),
    e3 = e3 + norm(vgg_get_nonhomg(y2(:,k))-y2b3(:,k))^2 + ...
        + norm(vgg_get_nonhomg(y1(:,k))-y1b3(:,k))^2;
end
e3 = sqrt(e3);

%ANSWER: e1 = 1.9029 e3 = 1.9032 Almost the same

%% Set another element then the last to 1

A02 = A(:,2:9);
b2 = A(:,1);

z2 = -(A02'*A02)^-1*A02'*b2;
z2 = [1 ; z2];


H2min = [z2(1:3) z2(4:6) z2(7:9)];

H2 = H2min;

y2b2 = vgg_get_nonhomg(H2*y1);
y1b2 = vgg_get_nonhomg(inv(H2)*y2);


figure(5);imagesc(img1);
figure(6);imagesc(img2);


figure(5);hold('on');plot(y1(1,:),y1(2,:),'go')
figure(6);hold('on');plot(y2(1,:),y2(2,:),'go')

figure(5);plot(y1b3(1,:),y1b3(2,:),'rx');
figure(6);plot(y2b3(1,:),y2b3(2,:),'rx');


%This shows the same result as with 1 in the last element
% It doesn't matter which position is set to 1 but it's easier to do
% caluclations if the first/last is set to 1


%% 4.4 Hartley-transformation

diag(S)';
figure(7);plot(log(diag(S)),'o');
 
% ANSWER: Nja they do not really represent a well-defined solution

[y1t T1]= liu_preconditioning(y1);
[y2t T2]= liu_preconditioning(y2);

y1tMeanDistance = 0;

for i = 1:length(y1t)
    y1tMeanDistance = y1tMeanDistance + sqrt(y1t(1,i).^2+y1t(2,i).^2);
end



y1tMeanDistance = y1tMeanDistance/length(y1t);
    
% ANSWER: y1tMeanDistance = 1.4142 = sqrt(2) -> The condition is satisfied!
% mean(y1t(1,:)) ~ 0

%%

At = [];
N = 8;


for cnt=1:N
    u = [y1t(1,cnt) y2t(1,cnt)];
    v = [y1t(2,cnt) y2t(2,cnt)];
    
    At = [At; [u(1) 0 -u(1)*u(2) v(1) 0 -v(1)*u(2) 1 0 -u(2)];
            [0 u(1) -u(1)*v(2) 0 v(1) -v(1)*v(2) 0 1 -v(2)]];
end
        
[Ut St Vt]=svd(At);
Ht = reshape(Vt(:,end),3,3);

H4 = T1^-1*Ht*T2;
y2b4 = vgg_get_nonhomg(H4*y1);
y1b4 = vgg_get_nonhomg(inv(H4)*y2);


%% Print and plot At

diag(St)';
figure(5);plot(log(diag(St)),'o');

% ANSWER: YES!! Awesome values!

%%  Symmetric geometric error of H4

e4 = 0;
for k = 1:length(y1),
    e4 = e4 + norm(vgg_get_nonhomg(y2(:,k))-y2b4(:,k))^2 + ...
        + norm(vgg_get_nonhomg(y1(:,k))-y1b4(:,k))^2;
end
e4 = sqrt(e4);


%ANSWER: e4 = 654.2985 VERY big!! Is this correct??????????????

%% 4.5 Ground Truth

load -ascii H1to2p
y2bGT = vgg_get_nonhomg(H1to2p*y1);
y1bGT = vgg_get_nonhomg(inv(H1to2p)*y2);

eGTtoH1 = 0;

for k = 1:length(y1b),
    eGTtoH1 = eGTtoH1 + norm(y2b(:,k)-y2bGT(:,k))^2 + ...
        + norm(y1b(:,k)-y1bGT(:,k))^2;
end
eGTtoH1 = sqrt(eGTtoH1);


%ANSWER: eGTtoH1 = 4.8963

eGTtoH2 = 0;

for k = 1:length(y1b),
    eGTtoH2 = eGTtoH2 + norm(y2b2(:,k)-y2bGT(:,k))^2 + ...
        + norm(y1b2(:,k)-y1bGT(:,k))^2;
end
eGTtoH2 = sqrt(eGTtoH2);
    

eGTtoH4 = 0;

for k = 1:length(y1b),
    eGTtoH4 = eGTtoH4 + norm(y2b4(:,k)-y2bGT(:,k))^2 + ...
        + norm(y1b4(:,k)-y1bGT(:,k))^2;
end
eGTtoH4 = sqrt(eGTtoH4);

% ANSWER: eGTtoH2 = 4.8961, eGTtoH4 = 652.1485
% H2 is the closest to the ground truth

%% 4.6 Transformations of Lines

y1line = y1(1:3,1:2);

x1 = y1line(1:3,1);
x2 = y1line(1:3,2);

l1plu = x1*x2'-x2*x1';
l2plu = H1*l1plu*H1';

k1 = l1plu(2,3)/l1plu(1,3);
m1 = l1plu(1,2)/l1plu(1,3);
l1 = [k1; -1; m1];
figure(1);drawline(l1,'axis','xy');

k2 = l2plu(2,3)/l2plu(1,3);
m2 = l2plu(1,2)/l2plu(1,3)
l2 = [k2; -1; m2];

figure(2);drawline(l2,'axis','xy');

% SJUKT SNYGGT!!!! The lines transforms to the places where they should!

%% 

y1onLine = vgg_get_homg([33 279; 267 447; 6 259; 83 311; 193 395; 321 487; 371 522; 420 556]');
y2onLine = round(H1*y1onLine);

N = length(y1onLine);

ALine = [];

for cnt=1:N
    u = [y1onLine(1,cnt) y2onLine(1,cnt)];
    v = [y1onLine(2,cnt) y2onLine(2,cnt)];
    
    ALine = [ALine; [u(1) 0 -u(1)*u(2) v(1) 0 -v(1)*u(2) 1 0 -u(2)];
            [0 u(1) -u(1)*v(2) 0 v(1) -v(1)*v(2) 0 1 -v(2)]];
end

[Ul Sl Vl] = svd(ALine);
HLine = reshape(Vl(:,end),3,3);

zLine = [HLine(1:3),H3(4:6),H3(7:9)]';

minError = abs(ALine*z3); %Corresponding calc for the minimal cas

diag(Sl)';
figure(8);plot(log(diag(Sl)),'o');

%ANSWER: The SVD have two points thats low... Not good


[y1onLinet T1Line]= liu_preconditioning(y1onLine);
[y2onLinet T2Line]= liu_preconditioning(y2onLine);

y1onLinetMeanDistance = 0;

for i = 1:length(y1onLinet)
    y1onLinetMeanDistance = y1onLinetMeanDistance + sqrt(y1onLinet(1,i).^2+y1onLinet(2,i).^2);
end



y1onLinetMeanDistance = y1onLinetMeanDistance/length(y1onLinet);


AonLinet = [];
N = 8;


for cnt=1:N
    u = [y1onLinet(1,cnt) y2onLinet(1,cnt)];
    v = [y1onLinet(2,cnt) y2onLinet(2,cnt)];
    
    AonLinet = [AonLinet; [u(1) 0 -u(1)*u(2) v(1) 0 -v(1)*u(2) 1 0 -u(2)];
            [0 u(1) -u(1)*v(2) 0 v(1) -v(1)*v(2) 0 1 -v(2)]];
end
        
[UonLinet SonLinet VonLinet]=svd(AonLinet);
HonLinet = reshape(VonLinet(:,end),3,3);

HLINE = T1Line^-1*HLine*T2Line;
y2bonLine = vgg_get_nonhomg(HLINE*y1onLine);
y1bonLine = vgg_get_nonhomg(inv(HLINE)*y2onLine);

diag(SonLinet)';
figure(9);plot(log(diag(SonLinet)),'o');

%%

img2tonLine=image_resample(double(img1),HLINE,640,800); %640, 800 size of image1
figure(10);imagesc(uint8(img2tonLine))

%Just a black figure, something is probably wrong

