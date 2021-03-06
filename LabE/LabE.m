%% LabE Principal Component Analysis

cd /site/edu/bb/MultidimensionalSignalAnalysis/ComputerExercises/ExerciseE/
addpath /edu/annhj876/Skola/TSBB06/LabE


%%

load A;
size(A);


% Compute the principal components of the dataset.
% Approch 1
C = A*A';                         %Compute the correlation matrix C
[e l] = eig(C);                   %Compute EVD of C
[PM p] = sort(diag(l),'descend');  %Sort the eigenvalues: largest first
PC = e(:,p);                      %Sort the eigenvectors in the same way
% Approch 2
[PC S] = svd(A);                  %Compute SVD of A
PM = diag(S);                     %magnitudes are given by the singular values

figure(1);
subplot(2,1,1);plot(PM,'o');
subplot(2,1,2);plot(log(PM),'o');

%How many of the principal components M do you consider to be significant
%for representing the signal?  Why?
% ANSWER: 4. F�r att det �r fyra egenv�rden som har betydligt st�rre
% magnitud �n de andra.

%% V�ljer optimalt M

M = 5;
figure(2);plot(diff(PM((M + 1):end)),'o');

% ANSWER: Vi valde M = 5, det �r inte samma v�rde som f�rut men n�stan

%% 

figure(3);plot(PC(:,1)'*A,PC(:,2)'*A,'o');axis('equal');
% ANSWER: The variance along the first principal axis is larger then the
% variance along the second.

%% 

m=mean(A,2);
A0 = A - m*ones(1,size(A,2));
figure(4);
subplot(2,1,1);plot(PM,'o');
subplot(2,1,2);plot(log(PM),'o');
figure(5);
subplot(4,1,1);plot(PC(:,1)'*A0,PC(:,2)'*A0,'o');axis('equal')
subplot(4,1,2);plot(PC(:,2)'*A0,PC(:,3)'*A0,'o');axis('equal')
subplot(4,1,3);plot(PC(:,3)'*A0,PC(:,4)'*A0,'o');axis('equal')
subplot(4,1,4);plot(PC(:,4)'*A0,PC(:,5)'*A0,'o');axis('equal')

%% VET INTE VAD VI G�R H�R. HJ�LP!!
error1 = cumsum(PC(1:M));
figure(6);plot(error1,'o');

%% 

figure(7);plot3(PC(:,1)'*A0,PC(:,2)'*A0,PC(:,3)'*A0,'o');axis('equal');

% ANSWER: The signal distrubution is Gaussian

% QUESTION: Try to suggest some other strategy for determining the number 
% of significant principal components M.

figure(8);mesh(A0)
% ANSWER: Yes we can make the same assumtions as when using plot3 but its a
% bit difficult to understand this plot...

%% 3 Principal component analysis of images

im_boat = double(imread('boat.png')); % Choose you image here!
image_rep(im_boat, 10,8) % SVAR: 8

im_baboon = double(imread('baboon.png')); % Choose you image here!
%image_rep(im_baboon, 6,8) % SVAR: 5

im_cameraman = double(imread('cameraman.png')); % Choose you image here!
%image_rep(im_cameraman, 18,8) % SVAR: 15

% ANSWER: Bra: Ej variationer i bilden. D�ligt: Kanter och variationer


%%

im = double(imread('middlebury.png')); % Choose you image here!
M = 6;

N = 8;
A = im2col(im,[N N],'distinct');
size(A)

C = A*A';        
[e l] = eig(C);
[PM p] = sort(diag(l),'descend');
PC = e(:,p);

[PC S] = svd(A);
PM = diag(S);

figure(11);
subplot(2,1,1);plot(PM,'o');
subplot(2,1,2);plot(log(PM),'o');

figure(12);colormap('gray');
c=PC(:,1:M)'*A;    %Compute coordinates from blocks
Arec=PC(:,1:M)*c;  %Reconstruct blocks from coordinates
imrec=col2im(Arec,[N N],size(im),'distinct');  %Reshape into image
imagesc(imrec);axis('off');     %Display image




im = double(imread('boat.png')); % Choose you image here!

figure(12);colormap('gray');
A = im2col(im,[N N],'distinct');
c=PC(:,1:M)'*A;    %Compute coordinates from blocks
Arec=PC(:,1:M)*c;  %Reconstruct blocks from coordinates
imrec=col2im(Arec,[N N],size(im),'distinct');  %Reshape into image
imagesc(imrec);axis('off');     %Display image

% ANSWER: There is nott a signifficant different. The six first principal
% vectors are lways almost the same

%%

figure(14);colormap('gray');
for cnt=1:6,
subplot(2,3,cnt);h=mesh(reshape(PC(:,cnt),N,N));
set(h,'edgecolor','black');axis([1 N 1 N -0.5 0.5]);
title(sprintf('principal component %d',cnt));
end

% ANSWER: plan, linj�r, exponensiellt kan ses fr�n log(EGENV�RDENA)

%% 3.2 Changing the block size

im_boat = double(imread('boat.png')); % Choose you image here!
image_rep(im_boat, 20,16)

% ANSWER: lower N better im quality but more block -> more data
% ANSWER: 18 is almost equalt to 8. when the blocksize is 16x16.
% Larger pixel blocks -> fewer blocks higher order principal components needed for the
% same detail


%% 3.3  PCA of synthetic images

im = double(imread('ploop.png')); % Choose you image here!
M = 64;

N = 8;
A = im2col(im,[N N],'distinct');
size(A)

C = A*A';        
[e l] = eig(C);
[PM p] = sort(diag(l),'descend');
PC = e(:,p);

[PC S] = svd(A);
PM = diag(S);

figure(11);
subplot(2,1,1);plot(PM,'o');
subplot(2,1,2);plot(log(PM),'o');

figure(12);colormap('gray');
c=PC(:,1:M)'*A;    %Compute coordinates from blocks
Arec=PC(:,1:M)*c;  %Reconstruct blocks from coordinates
imrec=col2im(Arec,[N N],size(im),'distinct');  %Reshape into image
imagesc(imrec);axis('off');     %Display image


figure(14);colormap('gray');
for cnt=1:6,
subplot(2,3,cnt);h=mesh(reshape(PC(:,cnt),N,N));
set(h,'edgecolor','black');axis([1 N 1 N -0.5 0.5]);
title(sprintf('principal component %d',cnt));
end


% ANSER: Yes there is a significant differance in the principal components 
% and principal values between ploop and the natural images. 

% ANSWER: No, the recunstruction is not good at all. Princible components
% are not linear . Looka at figure 11, the plot don't go down to zero.

%% Random image

im = rand(512,512);

M = 8;

N = 8;
A = im2col(im,[N N],'distinct');
size(A)

C = A*A';        
[e l] = eig(C);
[PM p] = sort(diag(l),'descend');
PC = e(:,p);

[PC S] = svd(A);
PM = diag(S);

figure(11);
subplot(2,1,1);plot(PM,'o');
subplot(2,1,2);plot(log(PM),'o');

figure(12);colormap('gray');
c=PC(:,1:M)'*A;    %Compute coordinates from blocks
Arec=PC(:,1:M)*c;  %Reconstruct blocks from coordinates
imrec=col2im(Arec,[N N],size(im),'distinct');  %Reshape into image
imagesc(imrec);axis('off');     %Display image


figure(14);colormap('gray');
for cnt=1:6,
subplot(2,3,cnt);h=mesh(reshape(PC(:,cnt),N,N));
set(h,'edgecolor','black');axis([1 N 1 N -0.5 0.5]);
title(sprintf('principal component %d',cnt));
end


% ANSWER: Yes, there is no consistency, that's what we expected.


%%

im = double(imread('boat.png')); % Choose you image here!

figure(12);colormap('gray');
A = im2col(im,[N N],'distinct');
c=PC(:,1:M)'*A;    %Compute coordinates from blocks
Arec=PC(:,1:M)*c;  %Reconstruct blocks from coordinates
imrec=col2im(Arec,[N N],size(im),'distinct');  %Reshape into image
imagesc(imrec);axis('off');     %Display image

% ANSWER: We can se what picture it is, but compared to how many M were
% needed previasly its bad. 

%% Use nataural image to reconstruct rand image

im = double(imread('middlebury.png')); % Choose you image here!
M = 8;

N = 8;
A = im2col(im,[N N],'distinct');
size(A)

C = A*A';        
[e l] = eig(C);
[PM p] = sort(diag(l),'descend');
PC = e(:,p);

[PC S] = svd(A);
PM = diag(S);

figure(11);
subplot(2,1,1);plot(PM,'o');
subplot(2,1,2);plot(log(PM),'o');

figure(12);colormap('gray');
c=PC(:,1:M)'*A;    %Compute coordinates from blocks
Arec=PC(:,1:M)*c;  %Reconstruct blocks from coordinates
imrec=col2im(Arec,[N N],size(im),'distinct');  %Reshape into image
imagesc(imrec);axis('off');     %Display image



im = rand(512,512);

figure(12);colormap('gray');
A = im2col(im,[N N],'distinct');
c=PC(:,1:M)'*A;    %Compute coordinates from blocks
Arec=PC(:,1:M)*c;  %Reconstruct blocks from coordinates
imrec=col2im(Arec,[N N],size(im),'distinct');  %Reshape into image
imagesc(imrec);axis('off');     %Display image

% ANSWER: The previous image were more detailed.