%MAPFoSt test using simulated images
raw=load('D:\Academics\Research\Seung Research\test images\out(0,0).mat');
raw=double(raw.out);
%raw=(raw(:,1:1024)+raw(:,1025:2048))/2;
raw=uint8(raw(:,1:1024));


psf=zeros(1024,1024);
psf(500:525,500:525)=1;
psf=psf/sum(sum(psf));
mtf=fft2(psf);

tic;
out=[];
for i=1:40
I1=raw;
I2=raw;
%test initial aberration
A=[i];

%% Initialize/calculate constants
%Calculate fft of two images
fI1=fft2(I1); %image should have dimension 2^n for faster FFT
fI2=fft2(I2);
height=size(I1,1);
width=size(I1,2);
FOV=8.511;
Acc=5;
PixSize = FOV/height; % um per pixel
A_max=80; %set max defocus and astigmatism to 80um-based of paper, needs to be changed
NA= 0.752 / (PixSize* (Acc*1000)^0.5);
%sigma=1; %estimated Gaussian noise (approximation for shot noise), rad/um (maybe calculate later)
sigma =mean([std(double(I1(:))), std(double(I2(:)))]); %sigma for real space
cutoffx=int32(floor(0.125*width)); %cutoff for k's used based on 25% k_nyquist, cycles/pixel
cutoffy=int32(floor(0.125*height)); %use for selecting subset of I/K e.g. K([1:cutoffy end+1-cutoffy:end],[1:cutoffx end+1-cutoffx:end])
%[Kx, Ky]=meshgrid((mod(0.5+[0:width-1]/width,1)-0.5)*(6.28/FOV),(mod(0.5+[0:height-1]/height,1)-0.5)*(6.28/FOV)); %units are rad/um?
[Kx, Ky]=meshgrid((circshift([0:width-1]/width,width/2,2)-0.5)*(6.28/FOV),(circshift([0:width-1]/width,width/2,2)-0.5)*(6.28/FOV)); %units are rad/um?

%[Kx,Ky] = meshgrid ([1:width]/(PixSize*width*1e-6),[1:height]/(PixSize*height*1e-6) );
%Kx=[1:width]*(6.28/FOV);
%Ky=[1:height]*(6.28/FOV);
%MTF=@(Kx,Ky,A) exp(-0.125*(NA^2)*(2*A(2)*(Kx.^2 - Ky.^2)*A(1) - A(3)*Kx.*Ky*A(1) + (Kx.^2+Ky.^2)*(A(2)^2 + A(3)^2 + A(1)^2)));

 MTF=@(Kx,Ky,A) exp(-0.125*(NA^2)*(Kx.^2+Ky.^2)*A^2);

% hardcoded test aberrations (defocus only)
T1=[15 ]; %defocus in [um]
T2=[-15 ];

noise=fft2(normrnd(0,sigma,height,width))+noise;
fI1=fI1.*MTF(Kx,Ky,A+T1);
%fI1=fI1.*mtf;
noise=fft2(normrnd(0,sigma,height,width));
fI2=fI2.*MTF(Kx,Ky,A+T2);
I1=(ifft2(fI1));
%I1=uint8(255*I1/max(max(I1)));
I2=(ifft2(fI2));
%I2=uint8(255*I2/max(max(I2)));

% close all
% figure;
% imshow([fftshift(fI1) fftshift(fI2)]);
% figure;
% imshow([I1 I2]);
% 
% fOmap=(fI1.*MTF(Kx,Ky,A+T1)+fI2.*MTF(Kx,Ky,A+T2))./(MTF(Kx,Ky,A+T1).^2+MTF(Kx,Ky,A+T2).^2);
% Omap=ifft2(fOmap);
% Omap=uint8(255*Omap/max(max(Omap)));
% figure;
% imshow(Omap);
% 
% MAP(A,I1,I2,T1,T2,FOV,Acc);
% checkgrad('MAP', randn(1,1), 1e-5,I1,I2,T1',T2',FOV,Acc);

O=minimize([0],@MAP,5,I1,I2,T1,T2,FOV,Acc);
out=[out O];
end
toc;
plot(1:20,out,1:20,1:20,':');
% set new WD/Stig from algorithm
% z=O(1);
% a_on=O(2);
% a_diag=O(3);