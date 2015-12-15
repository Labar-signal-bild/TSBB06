%% Lab F TSBB06 Over-sampling and the Discrete Wavelet Transform

cd /site/edu/bb/MultidimensionalSignalAnalysis/ComputerExercises/ExerciseF/
addpath /edu/annhj876/Skola/TSBB06/LabF

%% 3.1 Down-sampling and reconstruction of a discrete signal

fprintf('\n\n ------- 3.1 Down-sampling and reconstruction of a discrete signal -------\n\n')

s=randn(1024,1);
S=fftshift(fft(s));
u=(-512:511)'*pi/512;
Rect8=(abs(u)<pi/8);
Sbl=S.*Rect8; %Create a pi/4-band-limited transform
sbl=ifft(ifftshift(Sbl));
figure(1);subplot(2,1,1);plot(0:1023,sbl);
title('band-limited signal');
subplot(2,1,2);plot(u,abs(Sbl));
title('Fourier transform of band-limited signal');


downsampl8=sbl(1:8:end);

% Insert 7 zeros inbetween the downsampled signal to make it of the same
% length to make it more comparable
upsampl8=zeros(1024,1);
upsampl8(1:8:end)=downsampl8;
figure(2);subplot(2,1,1);plot(0:1023,upsampl8);
title('The signal down-sampled with a factor 8');
subplot(2,1,2);plot(u,abs(fftshift(fft(upsampl8))));
title('Fourier transform of down-sampled signal');

% Construct reconstuction filters
rect8=ifft(fftshift(8*Rect8));
figure(3);
subplot(2,1,1);plot(-512:511,ifftshift(rect8));
title('Reconstructing filter rect8');
subplot(2,1,2);plot(u,abs(fftshift(fft(rect8))));

for ix=0:127,
    B8(:,ix+1)=circshift(rect8,8*ix);
end
figure(4);plot(0:1023,B8(:,:));

% ANSWER: YES, the reconstructing filters are in accordance with what we
% thought.

% Reconstruct sbl
sblrec8=B8*upsampl8(1:8:end);
figure(5);
subplot(2,1,1);plot(0:1023,sblrec8);
title('Reconstructed signal from factor 8 down-sampled signal');
subplot(2,1,2);plot(u,abs(fftshift(fft(sblrec8))));
title('Fourier transform of reconstructed signal');
fprintf('Reconstruction from factor 8 down-sampled signal (no noise)\n\n');
fprintf('Reconstruction error: %e\n',norm(sbl-sblrec8));
fprintf('Reconstruction noise mean: %e\n',mean(sbl-sblrec8));
fprintf('Reconstruction noise std: %e\n\n',std(sbl-sblrec8));

%% 3.2 Sampling noise
fprintf('\n\n ------- 3.2 Sampling noise -------\n\n')

sigma=0.1;
downsampl8n=sbl(1:8:end)+sigma*randn(128,1);
upsampl8n=zeros(1,1024);
upsampl8n(1:8:end)=downsampl8n;
figure(6);
subplot(2,1,1);plot(0:1023,upsampl8n);
title('The signal down-sampled with a factor 8 with noise');
subplot(2,1,2);plot(u,abs(fftshift(fft(upsampl8n))));
title('Fourier transform of down-sampled signal with noise');

sblrec8n=B8*downsampl8n;
figure(7);subplot(2,1,1);plot(0:1023,sblrec8n);
title('Reconstructed signal from factor 8 down-sampled signal with noise');
subplot(2,1,2);plot(u,abs(fftshift(fft(sblrec8n))));
title('Fourier transform of reconstructed signal with noise');
err=sbl-sblrec8n;
fprintf('Reconstruction from factor 8 down-sampled signal with noise\n');
fprintf('Reconstruction noise error: mean %f std %f\n',mean(err),std(err));

% ANSWER: Mean~0 STD~0.1 same as prep exersice 5

%% 3.3    Over-sampling, reconstruction with a basis
fprintf('\n\n ------- 3.3 Over-sampling, reconstruction with a basis -------\n\n')

downsampl4n=sbl(1:4:end)+sigma*randn(256,1);
upsampl4n=zeros(1,1024);
upsampl4n(1:4:end)=downsampl4n;
figure(8);
subplot(2,1,1);plot(0:1023,upsampl4n);
title('The signal down-sampled with a factor 4 with noise');
subplot(2,1,2);plot(u,abs(fftshift(fft(upsampl4n))));
title('Fourier transform of down-sampled signal with noise');

Rect4=(abs(u)<=pi/4);
rect4=ifft(fftshift(4*Rect4));
figure(9);
subplot(2,1,1);plot(-512:511,ifftshift(rect4));
title('Reconstructing filter rect4');
subplot(2,1,2);plot(u,abs(fftshift(fft(rect4))));

for ix=0:255,
    B4(:,ix+1)=circshift(rect4,4*ix);
end
figure(10);plot(0:1023,B4(:,[1 2 3])');


% ANSWER: Same as before but less space between them.

sblrec4b=B4*downsampl4n;
figure(11);subplot(2,1,1);plot(0:1023,sblrec4b);
title('Reconstructed signal from factor 4 down-sampled signal (basis)');
subplot(2,1,2);plot(u,abs(fftshift(fft(sblrec4b))));
title('Fourier transform of reconstructed signal (basis)');
err=sbl-sblrec4b;
fprintf('Reconstruction from factor 4 down-sampled signal (basis)\n');
fprintf('Reconstruction noise error: mean %f std %f\n',mean(err),std(err));

% ANSWER: Mean~0 std~0.1

% Check if B4 is orthogonal
[U S V] = svd(B4);
diagS = diag(S);
zeros_S = find(diagS == 0); % zeros_S empty -> B4 orthgonal


%% 3.4 Over-sampling, reconstruction with a frame
fprintf('\n\n ------- 3.4 Over-sampling, reconstruction with a frame -------\n\n')

for ix=0:255,
    F8(:,ix+1)=circshift(rect8,4*ix);
end
figure(13);plot(0:1023,F8(:,[1 2 3])');

% ANSWER: The functions are more curved then before. The phase shift is not
% half a phase so the functions are not completely out of phase.


sblrec4f=0.5*F8*downsampl4n;
figure(15);
subplot(2,1,1);plot(0:1023,sblrec4f);
title('Reconstructed signal from factor 4 down-sampled signal (frame)');
subplot(2,1,2);plot(u,abs(fftshift(fft(sblrec4f))));
title('Fourier transform of reconstructed signal (frame)');
err=sbl-sblrec4f;
fprintf('Reconstruction from factor 4 down-sampled signal (frame)\n');
fprintf('Reconstruction noise error: mean %f std %f\n',mean(err),std(err));

% ANSWER: We get less noise since we don't use all the noise as we add to
% the signal.

% Answer preparatory 9: By taking the svd and checking so that 
% the diagonal is non zero.

% Check frame F8
[U S V] = svd(F8);
diagS  = diag(S);
zeros_S = find( diagS == 0); % zeros_S is empty -> F8 correct frame


%% 
% -------------------------------------------------------------------
% -------------------------------------------------------------------
% -------------------------------------------------------------------
% -------------------------------------------------------------------

%% 4.1 The two-channel filter bank
fprintf('\n\n ------- 4.1 The two-channel filter bank -------\n\n')

% Create random signal
l=256;
u=((-l/2):(l/2-1))*2*pi/l;
s0=rand(1,l);
S0=fftshift(fft(ifftshift(s0)));
S=S0.*(abs(u)<pi/4);
s=real(ifftshift(ifft(fftshift(S))));
figure(1);plot(0:255,s);
title('Input signal');

% Create filters from the db3 family
[h0 g0 h1 g1]=wfilters('db3')

% Look at the Fourier transforms of the filters
flength=length(h0);
u1=((-flength/2):(flength/2-1))*2*pi/flength;
figure(2);
subplot(4,1,1);plot(u1,abs(fftshift(fft(h0))));title('FT of h0');
subplot(4,1,2);plot(u1,abs(fftshift(fft(g0))));title('FT of g0');
subplot(4,1,3);plot(u1,abs((fftshift(fft(h1)))));title('FT of h1');
subplot(4,1,4);plot(u1,abs(fftshift(fft(g1))));title('FT of g1');

% ANSWER: h0,h1 are low-pass filters. g0,g1 are high-pass filters.
% ANSWER: Yes the coefficiants are as expected. fft(h0(end:-1:1)) =
% fft(h1), coefficienterna i h1 och  h0 är varandras spegling


dwtmode('per');  %Set periodic mode of filtering operations
[a d]=dwt(s,h0,g0);
figure(3);
subplot(2,1,1);plot(a);title('a');
subplot(2,1,2);plot(d);title('d');

% ANSWER: a är en approximation av signalen, de likaner varandratydligt. d
% är detaljerna, dess amplitud är mycket mindre än för signalen.

% Reconstruct
srec1=idwt(a,d,h1,g1);
figure(4);plot(0:255,srec1);
title('Reconstructed signal');

% Check if s and srec1 are equal
error = sum((s-srec1).^2); % Almost zero

%% 4.2 Multi-level filter bank
fprintf('\n\n ------- 4.2 Multi-level filter bank -------\n\n')

% Apply the filter on a0 instead
[a1 d1]=dwt(s,h0,g0);
[a2 d2]=dwt(a1,h0,g0);
figure(5);
subplot(3,1,1);plot(a2);title('a2');
subplot(3,1,2);plot(d2);title('d2');
subplot(3,1,3);plot(d1);title('d1');


a1rec=idwt(a2,d2,h1,g1);
srec2=idwt(a1rec,d1,h1,g1);
figure(6);plot(srec2);title('Reconstructed signal');

error = sum((s-srec2).^2); % A bit bigger then before but still really good

%% Iterations

% Code snippet (A)
N=5;ad=s;p=length(s);figure(7);
for cnt=1:N,
[a d]=dwt(ad(1:p),h0,g0);
ad(1:p)=[a d];

subplot(N+1,1,N+2-cnt);
plot(d);title(sprintf('details level %d',cnt));
p=p/2;
end
subplot(N+1,1,1);plot(a);
title(sprintf('approximation level %d',cnt));
figure(8);plot(ad);
title('Concatenated approximation and details');


% Reconstruct
% Code snippet (B)
for cnt=1:N,
ad(1:(2*p))=idwt(ad(1:p),ad((p+1):(2*p)),h1,g1);
p=2*p;
end
figure(9);plot(ad);title('Reconstructed signal');

error = sum((s-ad).^2); % Sämre men fortfarande väldigt bra

%% 4.3 Simple signal compression
fprintf('\n\n ------- 4.3 Simple signal compression -------\n\n')

% Code snippet (A)
N=3;ad=s;p=length(s);figure(7);
for cnt=1:N,
[a d]=dwt(ad(1:p),h0,g0);
ad(1:p)=[a d];

subplot(N+1,1,N+2-cnt);
plot(d);title(sprintf('details level %d',cnt));
p=p/2;
end
subplot(N+1,1,1);plot(a);
title(sprintf('approximation level %d',cnt));
figure(10);plot(ad);
title('Concatenated approximation and details');

p_save = p;
q = zeros((N+1),2);
q(:,1) = 4;
ad_test = ad;
ad_save = ad;

q(1,2) = max(abs(ad_test(1:p)));
ad_test = ad_test(p+1:end);

for i = 1:N
   q(i+1,2) = max(abs(ad_test(1:p)));
   ad_test = ad_test(p+1:end);
   p = 2*p; 
end

[ad bps] = quantisead(ad,q); %Replace the channels with quantised values

%%

p= p_save;
% Reconstruct
% Code snippet (B)
for cnt=1:N,
ad(1:(2*p))=idwt(ad(1:p),ad((p+1):(2*p)),h1,g1);
p=2*p;
end

error = sum((s-ad).^2);

figure(100);
subplot(3,1,1);plot(s);title('Original signal');
subplot(3,1,2);plot(ad);title('Reconstructed signal');
subplot(3,1,3);plot(s-ad);title('Difference');
SNR=log10(max(s)/std(ad-s))*20
bps

%% 4.3 My Quantisizer 
fprintf('\n\n ------- 4.3 My Quantisizer -------\n\n')

l=256;
u=((-l/2):(l/2-1))*2*pi/l;
s0=rand(1,l);
S0=fftshift(fft(ifftshift(s0)));
S=S0.*(abs(u)<pi/4);
s=real(ifftshift(ifft(fftshift(S))));

% Code snippet (A)
N=3;ad=s;p=length(s);figure(7);
for cnt=1:N,
[a d]=dwt(ad(1:p),h0,g0);
ad(1:p)=[a d];

subplot(N+1,1,N+2-cnt);
plot(d);title(sprintf('details level %d',cnt));
p=p/2;
end
subplot(N+1,1,1);plot(a);
title(sprintf('approximation level %d',cnt));
figure(10);plot(ad);
title('Concatenated approximation and details');

p_save = p;
q = zeros((N+1),2);

ad_test = ad;
ad_save = ad;
energy_vec = [];

q(1,2) = max(abs(ad_test(1:p)));
energy_vec = sum(ad_test(1:p))/q(1,2);
ad_test = ad_test(p+1:end);


for i = 1:N
   q(i+1,2) = max(abs(ad_test(1:p)));
   energy_vec = [energy_vec sum(ad_test(1:p))/q(i+1,2)];
   ad_test = ad_test(p+1:end);
   p = 2*p; 
end


energy_vec;
% q(:,1) = make_q(energy_vec,4);
q(:,1) = [6 3 2 1];
%q(:,1) = [6 3  1 1];

[ad bps] = quantisead(ad,q); %Replace the channels with quantised values



p= p_save;
% Reconstruct
% Code snippet (B)
for cnt=1:N,
ad(1:(2*p))=idwt(ad(1:p),ad((p+1):(2*p)),h1,g1);
p=2*p;
end



SNR=log10(max(s)/std(ad-s))*20
bps
