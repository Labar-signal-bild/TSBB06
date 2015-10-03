%% Lab1 uppgift 4

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

%A = A(1:8,:);

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

%ANSWER: Very close, about half a pixel wrong.
%ANSWER: Yes, A*z ~ 0.1 - 0.8


%%  Symmetric geometric error

e1 = 0;
for k = 1:length(y1),
    e1 = e1 + norm(vgg_get_nonhomg(y2(:,k))-y2b(:,k))^2 + ...
        + norm(vgg_get_nonhomg(y1(:,k))-y1b(:,k))^2;
end
e1 = sqrt(e1);

%ANSWER: e1 = 1.9029 Kind of the same as expected


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
%located where the theres a colorchange. In the bottom right corner0 of 
%figure 4 there's a lot of color, that's because in image1 this isn't even shown.

%There is a car in the picture who is moving. So in the second image the
%car is covering more of the picture then in image1


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
%Rational reason to set a particular element to 1??????????????????

%% 4.3 Homogeneous solution
%A = A(1:8,:); %A0
[U S V] = svd(A);
H3 = reshape(V(:,end),3,3);

z3 = [H3(1:3),H3(4:6),H3(7:9)]';

minError = abs(A*z3); %Corresponding calc for the minimal case

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

%% 4.4 Hartley-transformation


 










