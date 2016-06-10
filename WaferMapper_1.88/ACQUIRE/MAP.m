function [f, df]=MAP(A,fI1,fI2,T1,T2,FOV,height,width)
%MAP function to maximize with respect to aberration vector, A

%% Initialize/calculate constants
A_max=80; %set max defocus and astigmatism to 80um-based of paper, needs to be changed
NA=5; %set NA to 5 mrad-based of paper, needs to be changed
%NA= 0.752 / (PixSize* (Acc*1000)^0.5);
sigma=1; %estimated Gaussian noise (approximation for shot noise), rad/um (maybe calculate later)
cutoffx=int(floor(0.125*width)); %cutoff for k's used based on 25% k_nyquist, cycles/pixel
cutoffy=int(floor(0.125*height)); %use for selecting subset of I/K e.g. I([1:cutoffy end+1-cutoffy:end],[1:cutoffx end+1-cutoffx:end])
[Kx Ky]=meshgrid([1:width]*(6.28/FOV),[1:height]*(6.28/FOV));
%Kx=[1:width]*(6.28/FOV);%units are rad/um?
%Ky=[1:height]*(6.28/FOV);

%% Setup probability functions for MAPFoSt
p_A=@(A)  (max(abs(A(1)),abs(A(2)),abs(A(3)))>=A_max)/((2*A_max)^3); %uniform distribtion over cube of side 2*A_max
pA=p_A(A);
%MTF=@(kx,ky,A) exp(-0.125*NA^2*(2*A(2)*(kx^2-ky^2)*A(1)-A(3)*kx*ky*A(1)+(kx^2+ky^2)*(A(1)^2+A(2)^2+A(3)^2))); %check for typos
MTF=@(Kx,Ky,A) exp(-0.125*(NA^2)*(2*A(2)*(Kx.^2 - Ky.^2)*A(1) - A(3)*Kx.*Ky*A(1) + (Kx.^2+Ky.^2).*(A(2)^2 + A(3)^2 + A(1)^2)));
f=-log(p_A(A))+sum(sum( ((I2.*MTF(Kx,Ky,A.+T1) - I1.*MTF(Kx,Ky,A.+T2)).^2) ./ (2*sigma^2*(MTF(Kx,Ky,A.+T1).^2+MTF(Kx,Ky,A.+T2).^2));

%MTF equation used to calculate partial derivatives in Wolfram
%exp(-0.125*(N^2)*(2a*(x^2-y^2)*z-b*x*y*z+(x^2+y^2)(a^2+b^2+z^2)))
const=-0.125*(NA^2)*exp(-0.125*(NA^2)*(2*A(2)*(Kx.^2 - Ky.^2)*A(1) - A(3)*Kx.*Ky*A(1) + (Kx.^2+Ky.^2).*(A(2)^2 + A(3)^2 + A(1)^2)));
dMTFd1=@(Kx,Ky,A) (2*A(2)*(Kx.^2 - Ky.^2) + 2*A(1)*(Kx.^2 + Ky.^2) - A(3)*Kx.*Ky).*const;
dMTFd2=@(Kx,Ky,A) (2*A(2)*(Kx.^2 + Ky.^2) + 2*A(1)*(Kx.^2 - Ky.^2)).*const;
dMTFd3=@(Kx,Ky,A) (2*A(3)*(Kx.^2 + Ky.^2) - A(1)*Kx.*Ky).*const;

Hi=((I2.*MTF(Kx,Ky,A.+T1) - I1.*MTF(Kx,Ky,A.+T2)).^2;
Lo=(2*sigma^2*(MTF(Kx,Ky,A.+T1).^2+MTF(Kx,Ky,A.+T2).^2));

df(1)= sum(sum( (Lo.*(I2.*dMTFd1(Kx,Ky,A.+T1) - I1.*dMTFd1(Kx,Ky,A.+T2)) - Hi.*(4*sigma^2*(MTF(Kx,Ky,A.+T1).*dMTFd1(Kx,Ky,A.+T1).+MTF(Kx,Ky,A.+T2).*dMTFd1(Kx,Ky,A.+T2)))) / Lo.^2 ));
df(2)= sum(sum( (Lo.*(I2.*dMTFd2(Kx,Ky,A.+T1) - I1.*dMTFd2(Kx,Ky,A.+T2)) - Hi.*(4*sigma^2*(MTF(Kx,Ky,A.+T1).*dMTFd2(Kx,Ky,A.+T1).+MTF(Kx,Ky,A.+T2).*dMTFd2(Kx,Ky,A.+T2)))) / Lo.^2 ));
df(3)= sum(sum( (Lo.*(I2.*dMTFd3(Kx,Ky,A.+T1) - I1.*dMTFd3(Kx,Ky,A.+T2)) - Hi.*(4*sigma^2*(MTF(Kx,Ky,A.+T1).*dMTFd3(Kx,Ky,A.+T1).+MTF(Kx,Ky,A.+T2).*dMTFd3(Kx,Ky,A.+T2)))) / Lo.^2 ));

end