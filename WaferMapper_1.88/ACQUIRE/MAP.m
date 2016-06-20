function [f, df]=MAP(A,I1,I2,T1,T2,FOV,Acc)
%MAP function to maximize with respect to aberration vector, A
%T1 and T2 are the known test aberrations

%% Initialize/calculate constants
%Calculate fft of two images
fI1=fft2(I1); %image should have dimension 2^n for faster FFT
fI2=fft2(I2);
height=size(I1,1);
width=size(I1,2);
PixSize = FOV/height; % um per pixel
A_max=80; %set max defocus and astigmatism to 80um-based of paper, needs to be changed
%NA=5; %set NA to 5 mrad-based of paper, needs to be changed
NA= 0.752 / (PixSize* (Acc*1000)^0.5); %confirm units for PixSize
%sigma=1; %estimated Gaussian noise (approximation for shot noise), rad/um (maybe calculate later)
sigma =mean([std(double(I1(:))), std(double(I2(:)))]);
%cutoffx=int32(floor(0.125*width)); %cutoff for k's used based on 25% k_nyquist, cycles/pixel
%cutoffy=int32(floor(0.125*height)); %use for selecting subset of I/K e.g. K([1:cutoffy end+1-cutoffy:end],[1:cutoffx end+1-cutoffx:end])
[Kx, Ky]=meshgrid((circshift([0:width-1]/width,width/2,2)-0.5)*(6.28/FOV),(circshift([0:width-1]/width,width/2,2)-0.5)*(6.28/FOV)); %units are rad/um?
%Kx=[1:width]*(6.28/FOV);
%Ky=[1:height]*(6.28/FOV);

%% Setup probability functions for MAPFoSt
p_A=@(A)  (max([abs(A(1)),abs(A(2)),abs(A(3))])<=A_max)/((2*A_max)^3); %uniform distribtion over cube of side 2*A_max
%MTF=@(kx,ky,A) exp(-0.125*NA^2*(2*A(2)*(kx^2-ky^2)*A(1)-A(3)*kx*ky*A(1)+(kx^2+ky^2)*(A(1)^2+A(2)^2+A(3)^2))); %check for typos
MTF=@(Kx,Ky,A) exp(-0.125*(NA^2)*(2*A(2)*(Kx.^2 - Ky.^2)*A(1) - A(3)*Kx.*Ky*A(1) + (Kx.^2+Ky.^2).*(A(2)^2 + A(3)^2 + A(1)^2)));
f=-log(p_A(A))+sum(sum( (abs(fI2.*MTF(Kx,Ky,A+T1) - fI1.*MTF(Kx,Ky,A+T2)).^2) ./ (2*sigma^2*(MTF(Kx,Ky,A+T1).^2+MTF(Kx,Ky,A+T2).^2)+1e-10) ));

%MTF equation used to calculate partial derivatives in Wolfram
%exp(-0.125*(N^2)*(2a*(x^2-y^2)*z-b*x*y*z+(x^2+y^2)(a^2+b^2+z^2)))
const=-0.125*(NA^2)*exp(-0.125*(NA^2)*(2*A(2)*(Kx.^2 - Ky.^2)*A(1) - A(3)*Kx.*Ky*A(1) + (Kx.^2+Ky.^2).*(A(2)^2 + A(3)^2 + A(1)^2)));
dMTFd1=@(Kx,Ky,A) (2*A(2)*(Kx.^2 - Ky.^2) + 2*A(1)*(Kx.^2 + Ky.^2) - A(3)*Kx.*Ky).*const;
dMTFd2=@(Kx,Ky,A) (2*A(2)*(Kx.^2 + Ky.^2) + 2*A(1)*(Kx.^2 - Ky.^2)).*const;
dMTFd3=@(Kx,Ky,A) (2*A(3)*(Kx.^2 + Ky.^2) - A(1)*Kx.*Ky).*const;

Hi=((fI2.*MTF(Kx,Ky,A+T1) - fI1.*MTF(Kx,Ky,A+T2)).^2);
Lo=(2*sigma^2*(MTF(Kx,Ky,A+T1).^2+MTF(Kx,Ky,A+T2).^2)+1e-10);

df(1)= sum(sum( (Lo.*(fI2.*dMTFd1(Kx,Ky,A+T1) - fI1.*dMTFd1(Kx,Ky,A+T2)) - Hi.*(4*sigma^2*(MTF(Kx,Ky,A+T1).*dMTFd1(Kx,Ky,A+T1)+MTF(Kx,Ky,A+T2).*dMTFd1(Kx,Ky,A+T2)))) / Lo.^2 ));
df(2)= sum(sum( (Lo.*(fI2.*dMTFd2(Kx,Ky,A+T1) - fI1.*dMTFd2(Kx,Ky,A+T2)) - Hi.*(4*sigma^2*(MTF(Kx,Ky,A+T1).*dMTFd2(Kx,Ky,A+T1)+MTF(Kx,Ky,A+T2).*dMTFd2(Kx,Ky,A+T2)))) / Lo.^2 ));
df(3)= sum(sum( (Lo.*(fI2.*dMTFd3(Kx,Ky,A+T1) - fI1.*dMTFd3(Kx,Ky,A+T2)) - Hi.*(4*sigma^2*(MTF(Kx,Ky,A+T1).*dMTFd3(Kx,Ky,A+T1)+MTF(Kx,Ky,A+T2).*dMTFd3(Kx,Ky,A+T2)))) / Lo.^2 ));

end