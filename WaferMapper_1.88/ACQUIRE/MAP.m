function [f, df]=MAP(A,I1,I2,T1,T2,FOV,Acc,single)
%MAP function to maximize with respect to aberration vector, A
%T1 and T2 are the known test aberrations
gaussian =1;
%% Initialize/calculate constants
%Calculate fft of two images
fI1=fft2(I1); %image should have dimension 2^n for faster FFT
fI2=fft2(I2);
height=size(I1,1);
width=size(I1,2);
PixSize = FOV/height; % um per pixel
A_max=80; %set max defocus and astigmatism to 80um-based of paper, needs to be changed
A_sigma=80; %sigma for gaussian prior
%NA=5; %set NA to 5 mrad-based of paper, needs to be changed
NA= 0.752 / (PixSize* (Acc*1000)^0.5); %confirm units for PixSize
%sigma=10; %estimated Gaussian noise (approximation for shot noise), rad/um (maybe calculate later)
sigma =mean([std(double(I1(:))), std(double(I2(:)))]);
cutoffx=int32(floor(0.125*width)); %cutoff for k's used based on 25% k_nyquist, cycles/pixel
cutoffy=int32(floor(0.125*height)); %use for selecting subset of I/K e.g. K([1:cutoffy end+1-cutoffy:end],[1:cutoffx end+1-cutoffx:end])
[Kx, Ky]=meshgrid((mod(0.5+[0:width-1]/width,1)-0.5)*(6.28/FOV),(mod(0.5+[0:height-1]/height,1)-0.5)*(6.28/FOV));
% determine which wave vectors to use
Kx=Kx([1:cutoffy end+1-cutoffy:end],[1:cutoffx end+1-cutoffx:end]);
Ky=Ky([1:cutoffy end+1-cutoffy:end],[1:cutoffx end+1-cutoffx:end]);
fI1=fI1([1:cutoffy end+1-cutoffy:end],[1:cutoffx end+1-cutoffx:end]);
fI2=fI2([1:cutoffy end+1-cutoffy:end],[1:cutoffx end+1-cutoffx:end]);
% Kx=Kx(1,:);
% Ky=Ky(1,:);
% fI1=fI1(1,:);
% fI2=fI2(1,:);
%precompute some matrices
Kx2=Kx.^2;
Ky2=Ky.^2;

if single
    %for single aberration mode
    Kx2Ky2=Kx2+Ky2;
    p_A=@(A)  (max(abs(A))<=A_max)/(2*A_max); %uniform distribtion over cube of side 2*A_max
    MTF=@(Kx,Ky,A) exp(-5*(NA^2)*(Kx2Ky2)*A^2);
    dMTFdz=@(Kx,Ky,A) -10*(NA^2)*(Kx2Ky2)*A.*exp(-5*(NA^2)*(Kx2Ky2)*A^2);
    MTF1=MTF(Kx,Ky,A+T1);
    MTF2=MTF(Kx,Ky,A+T2);
    temp=((fI2.*MTF1 - fI1.*MTF2).^2) ./ (2*sigma^2*(MTF1.^2+MTF2.^2)+1e-20); %this is the matrix of log likelihoods (ie the exponentiated part of eq 27)
    %df=sum(sum( (2*(fI2.*MTF1-fI1.*MTF2).*(fI2.*dMTFdz(Kx,Ky,A+T1)-fI1.*dMTFdz(Kx,Ky,A+T2)))./(2*sigma^2*(MTF1.^2+MTF2.^2)+1e-20)-(2*(MTF1.*dMTFdz(Kx,Ky,A+T1)+MTF2.*dMTFdz(Kx,Ky,A+T2)).*(fI2.*MTF1-fI1.*MTF2).^2)./(2*sigma^2*(MTF1.^2+MTF2.^2).^2+1e-20) ));
    df=(2*(MTF1.*dMTFdz(Kx,Ky,A+T2)-MTF2.*dMTFdz(Kx,Ky,A+T1)).*(fI2.*MTF2+fI1.*MTF1).*(fI1.*MTF2-fI2.*MTF1))./(2*sigma^2*(MTF1.^2+MTF2.^2).^2+1e-20) ;
    df=sum(sum( df(temp>0) ));
    if gaussian
        p_A=@(A)  (1/(2.5066*A_sigma))*exp(-A^2/(2*A_sigma^2)); %gaussian prior
        df=df+A/A_sigma^2;
    end
    f=-log(p_A(A))+sum(sum( temp(temp>0) ));
else
    %% Setup probability functions for MAPFoSt
    p_A=@(A)  (max([abs(A(1)),abs(A(2)),abs(A(3))])<=A_max)/((2*A_max)^3); %uniform distribtion over cube of side 2*A_max
    MTF=@(Kx,Ky,A) exp(-0.125*(NA^2)*(2*A(2)*(Kx2 - Ky2)*A(1) - A(3)*Kx.*Ky*A(1) + (Kx2+Ky2)*(A(2)^2 + A(3)^2 + A(1)^2)));
    MTF1=MTF(Kx,Ky,A+T1);
    MTF2=MTF(Kx,Ky,A+T2);
    %MTF equation used to calculate partial derivatives in Wolfram
    %exp(-0.125*(N^2)*(2a*(x^2-y^2)*z-b*x*y*z+(x^2+y^2)(a^2+b^2+z^2)))
    dMTFd1=@(Kx,Ky,A) -0.125*(NA^2)*(2*A(2)*(Kx2 - Ky2) + 2*A(1)*(Kx2 + Ky2) - A(3)*Kx.*Ky).*exp(-0.125*(NA^2)*(2*A(2)*(Kx2 - Ky2)*A(1) - A(3)*Kx.*Ky*A(1) + (Kx2+Ky2)*(A(2)^2 + A(3)^2 + A(1)^2)));
    dMTFd2=@(Kx,Ky,A) -0.125*(NA^2)*(2*A(2)*(Kx2 + Ky2) + 2*A(1)*(Kx2 - Ky2)).*exp(-0.125*(NA^2)*(2*A(2)*(Kx2 - Ky2)*A(1) - A(3)*Kx.*Ky*A(1) + (Kx2+Ky2)*(A(2)^2 + A(3)^2 + A(1)^2)));
    dMTFd3=@(Kx,Ky,A) -0.125*(NA^2)*(2*A(3)*(Kx2 + Ky2) - A(1)*Kx.*Ky).*exp(-0.125*(NA^2)*(2*A(2)*(Kx.^2 - Ky.^2)*A(1) - A(3)*Kx.*Ky*A(1) + (Kx2+Ky2)*(A(2)^2 + A(3)^2 + A(1)^2)));

    % d/dx (a*f(x+i)-b*f(x+j))^2/(c*(f(x+i)^2+f(x+j)^2))
    % d/dx (((a+c*i)*f(x+g)-(b+d*i)*f(x+h))^2)/(s*(f(x+g)^2+f(x+h)^2))
    df(1)= sum(sum( (2*(MTF1.*dMTFd1(Kx,Ky,A+T2)-MTF2.*dMTFd1(Kx,Ky,A+T1)).*(fI2.*MTF2+fI1.*MTF1).*(fI1.*MTF2-fI2.*MTF1))./(2*sigma^2*(MTF1.^2+MTF2.^2).^2) ));
    df(2)= sum(sum( (2*(MTF1.*dMTFd2(Kx,Ky,A+T2)-MTF2.*dMTFd2(Kx,Ky,A+T1)).*(fI2.*MTF2+fI1.*MTF1).*(fI1.*MTF2-fI2.*MTF1))./(2*sigma^2*(MTF1.^2+MTF2.^2).^2) ));
    df(3)= sum(sum( (2*(MTF1.*dMTFd3(Kx,Ky,A+T2)-MTF2.*dMTFd3(Kx,Ky,A+T1)).*(fI2.*MTF2+fI1.*MTF1).*(fI1.*MTF2-fI2.*MTF1))./(2*sigma^2*(MTF1.^2+MTF2.^2).^2) ));
    f=-log(p_A(A))+sum(sum( ((fI2.*MTF1 - fI1.*MTF2).^2) ./ (2*sigma^2*(MTF1.^2+MTF2.^2)) ));
    %f=sum(sum( ((fI2.*MTF1 - fI1.*MTF2).^2) ./ (2*sigma^2*(MTF1.^2+MTF2.^2)) ));
end
end