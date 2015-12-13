%% TSBB06 Computer Exercise C Normalised Convolution

cd /site/edu/bb/MultidimensionalSignalAnalysis/ComputerExercises/ExerciseC/
addpath /edu/annhj876

%% Polynomial Expansion

% Signal s
k = 0:1:100;
s = sin(k/10);
figure(1);plot(s);

% Basis functions
x=(-3:3)';
b0=ones(7,1);
b1=x;
b2=x.^2;
figure(2);
subplot(4,1,1);plot(b0,'-o');
subplot(4,1,2);plot(b1,'-o');
subplot(4,1,3);plot(b2,'-o');

% Applicability function, use a gaussian that fits within 7 points
a = exp(-x.^2/4);
figure(2);
subplot(4,1,4);plot(a,'-o');

% Construct filters by multiplying the basis functions with the
% applicability function
f0 = b0.*a; f0 = f0(end:-1:1);
f1 = b1.*a; f1 = f1(end:-1:1);
f2 = b2.*a; f2 = f2(end:-1:1);
figure(3);
subplot(3,1,1);plot(f0,'-o');
subplot(3,1,2);plot(f1,'-o');
subplot(3,1,3);plot(f2,'-o');


% Convolve the signal with the filters
h0 = conv(s,f0,'same');
h1 = conv(s,f1,'same');
h2 = conv(s,f2,'same');

% Create the metric G
G0 = diag(a);
B = [b0 b1 b2];
G = B'*G0*B

% ANSWER: The basis are not orthonormal. When observing G we see that the 
% two diagonals are not zero, hence the basis are not all orthonormal to each other.  


% Compute proper coordinates
% The plot shows the elements of the coordinate vector c
c = inv(G)*[h0;h1;h2];
figure(4);
subplot(3,1,1);plot(c(1,:))
subplot(3,1,2);plot(c(2,:))
subplot(3,1,3);plot(c(3,:))

% Are the plots as expected? Hint:  polynomial expansion is related to a Taylor
% expansion  of  a  function,  as  described  in  preparatory  exercise  2.   
% In  what way?
% ANSWER: The first one looks like we expected, it should be the original
% signal which it is. The second should (the first order derivitive) be a
% scaled cosine, which it is. The last one a sine with amplitude 0.01


figure(5);
localsig=s(60-3:60+3);
reconsig=(B*c(:,60))';
diffsig=localsig-reconsig;
subplot(3,1,1);plot(localsig);
subplot(3,1,2);plot(reconsig);
subplot(3,1,3);plot(diffsig);

% The reconstructed signal is similar to the original signal. It's more
% simular in the middle because the applicability function is 1 in the
% middle and zero on the edges which means that it allows the signal to be
% not as simular on the edges. 

%% 4 Uncertain data

% his produces a signal vector svert where 20% of the samples are set to zero

cert = double(rand(1,101)>0.5);
scert = s.*cert;
figure(6);plot(scert);

% Convolve the signal with the filters
h02 = conv(scert,f0,'same');
h12 = conv(scert,f1,'same');
h22 = conv(scert,f2,'same');


% Compute proper coordinates
% The plot shows the elements of the coordinate vector c
c2 = inv(G)*[h02;h12;h22];
figure(7);
subplot(3,1,1);plot(c2(1,:))
subplot(3,1,2);plot(c2(2,:))
subplot(3,1,3);plot(c2(3,:))


figure(8);
localsig2=scert(60-3:60+3);
reconsig2=(B*c2(:,60))';
diffsig2=localsig2-reconsig2;
subplot(3,1,1);plot(localsig2);
subplot(3,1,2);plot(reconsig2);
subplot(3,1,3);plot(diffsig2);

% ANSWER: The result change alot when having an uncertainy of the signal.
% In the first case we can still see it looks like i sine but in the first
% derivitive it does not at all look the same as before, this since the
% signal change alot which make the derievitive change even more. 


%%

h0 = conv(scert,f0,'same'); 
h1 = conv(scert,f1,'same');

f11 = b0.*a.*b0; f11 = f11(end:-1:1);
f12 = b0.*a.*b1; f12 = f12(end:-1:1);
f22 = b1.*a.*b1; f22 = f22(end:-1:1);

G11 = conv(cert,f11,'same');
G12 = conv(cert,f12,'same');
G22 = conv(cert,f22,'same');

detG = G11.*G22-G12.^2;
c0 = (G22.*h0-G12.*h1)./detG; % inv(G)*h
c1 = (-G12.*h0+G11.*h1)./detG; 
figure(9);
subplot(2,1,1);plot(c0)
subplot(2,1,2);plot(c1)

% How do you interpret the computations that are made for c0
% and c1?
% ANSWER: It's the inverse of G times h that give the basis c0 c1

% Compare c0 and c1 from before: Now it has spikar which it didn't before. 


% Modify the signal by removing a higher rate of the samples
% (e.g. 50%) and see how this changes the result. How large rate of samples
% can be removed without too much distortion in the measured coordinates
% (c0, c1)?
% ANSWER: Ca 30% 

%% 5 Normalised averaging of images

im = double(imread('Scalespace0.png'));
figure(10);colormap(gray);imagesc(im);


cert = double(rand(size(im)) > 0.6); 
imcert = im.*cert;
figure(11);colormap(gray);imagesc(imcert);

% Applicability functian a is a low-pass filter
x = ones(7,1)*(-3:3)
y = x'
a = exp(-(x.^2+y.^2)/4);
figure(12);mesh(a);

% Low pass filtering the uncertian image with applicability function
imlp = conv2(imcert, a, 'same');
figure(13);colormap(gray);imagesc(imlp);

% ANSWER: Yes this did solve the problem with missing pixels since we
% smudge the image.


G = conv2(cert, a, 'same'); % We wrote this line
c = imlp./G;
figure(14);colormap(gray);imagesc(c);

% ANSWER: This reault is much better then before. G include information
% about where the missing pixels are and can from that reconstruct the
% image. Why does it not work with imcert??


% Change width of a: longer: More smudged
%                    shorter: Not all of the pixelvalues are getting a
%                    value


% ANSWER: When removing 97% of the pixels we can not longer see what the
% image show

% ANSWER: When 90% of the pixels are gone and a has the length of 7 there
% will be zeros in G which we then divide the imlp image with -> dividing
% by 0 give NaN result.

% ANSWER: When increasing the lenth of a we improve the image.