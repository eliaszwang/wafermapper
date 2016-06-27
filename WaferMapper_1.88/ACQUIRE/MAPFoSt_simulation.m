%MAPFoSt test using simulated images
close all
clear
raw=load('D:\Academics\Research\Seung Research\test images\out(0,0)');
raw=(raw.out);
raw1=(raw(:,1:1024));
raw2=(raw(:,1025:2048));
raw3=raw(1:512,1:512);
% psf=zeros(1024,1024);
% psf(500:525,500:525)=1;
% psf=psf/sum(sum(psf));
% mtf=fft2(psf);


single=1;
tic;
out=[];
r=1:10; %loop over different (actual) aberrations
for i=r
%% Initialize/calculate constants
%Calculate fft of two images
fI1=fft2((raw2)); %image should have dimension 2^n for faster FFT
fI2=fft2((raw2));
height=size(fI1,1);
width=size(fI1,2);
FOV=8.511;
Acc=5;
PixSize = FOV/height; % um per pixel
A_max=80; %set max defocus and astigmatism to 80um-based of paper, needs to be changed
NA= 0.752 / (PixSize* (Acc*1000)^0.5);
sigma=1; %estimated Gaussian noise (approximation for shot noise), rad/um (maybe calculate later)
%sigma =mean([std(double(I1(:))), std(double(I2(:)))]); %sigma for real space
cutoffx=int32(floor(0.125*width)); %cutoff for k's used based on 25% k_nyquist, cycles/pixel
cutoffy=int32(floor(0.125*height)); %use for selecting subset of I/K e.g. K([1:cutoffy end+1-cutoffy:end],[1:cutoffx end+1-cutoffx:end])
%[Kx, Ky]=meshgrid((mod(0.5+[0:width-1]/width,1)-0.5)*(6.28/FOV),(mod(0.5+[0:height-1]/height,1)-0.5)*(6.28/FOV)); %units are rad/um?
[Kx, Ky]=meshgrid((circshift([0:width-1]/width,width/2,2)-0.5)*(6.28/FOV),(circshift([0:width-1]/width,width/2,2)-0.5)*(6.28/FOV)); %units are rad/um?

if single
    %test initial aberration
    A=i;
    % hardcoded test aberrations (defocus only)
    T1=15; %defocus in [um]
    T2=-15;
    init=0;
    MTF=@(Kx,Ky,A) exp(-5*(NA^2)*(Kx.^2+Ky.^2)*A^2);
else
    %test initial aberration
    A=[i 3 5];
    % hardcoded test aberrations (defocus only)
    T1=[15 0 0]; %defocus in [um]
    T2=[-15 0 0];
    init=[0 0 0];
    MTF=@(Kx,Ky,A) exp(-0.125*(NA^2)*(2*A(2)*(Kx.^2 - Ky.^2)*A(1) - A(3)*Kx.*Ky*A(1) + (Kx.^2+Ky.^2)*(A(2)^2 + A(3)^2 + A(1)^2)));
end




noise=fft2(normrnd(0,sigma,height,width));
fI1=fI1.*MTF(Kx,Ky,A+T1)+noise;
noise=fft2(normrnd(0,sigma,height,width));
fI2=fI2.*MTF(Kx,Ky,A+T2)+noise;
I1=(ifft2(fI1));
%I1=uint8(255*I1/max(max(I1)));
I2=(ifft2(fI2));
%I2=uint8(255*I2/max(max(I2)));

% close all
% figure;
% imshow([fftshift(fI1) fftshift(fI2)]);
% figure;
% imshow([uint8(255*I1/max(max(I1))) uint8(255*I2/max(max(I2)))]);
% 
% fOmap=(fI1.*MTF(Kx,Ky,A+T1)+fI2.*MTF(Kx,Ky,A+T2))./(MTF(Kx,Ky,A+T1).^2+MTF(Kx,Ky,A+T2).^2);
% Omap=ifft2(fOmap);
% Omap=uint8(255*Omap/max(max(Omap)));
% figure;
% imshow(Omap);
% 
% MAP(A,I1,I2,T1,T2,FOV,Acc,single);
% checkgrad('MAP', randn(1,1), 1e-5,I1,I2,T1',T2',FOV,Acc,single);
p.length=10;
p.method='BFGS';
p.verbosity=0;
p.MFEPLS = 10;   % Max Func Evals Per Line Search
p.MSR = 100;                % Max Slope Ratio default
O=minimize(init,@MAP,p,I1,I2,T1,T2,FOV,Acc,single);
out=[out [A';O';MAP(A,I1,I2,T1,T2,FOV,Acc,single);MAP(O,I1,I2,T1,T2,FOV,Acc,single)]];
%plot 1-D MAP
fout=[];
dfout=[];
temp=-10:10;
for j=temp
    [f,df]=MAP(j,I1,I2,T1,T2,FOV,Acc,1);
    fout=[fout f ];
    dfout=[dfout df];
end
h=figure;
plot(temp,fout);
yyaxis right;
plot(temp,dfout);
title(['simulation -ln(P(A)) vs A for A=' num2str(A)]);
%saveas(h,['D:\Academics\Research\Seung Research\Analysis plots\simulation raw3 -ln(P(A)) vs A for A=' num2str(A) '.jpg']);

end
% MSE=immse(out(1,:),out(2,:));
% disp(['MSE: ' num2str(MSE)]);
toc;

if single
    figure;
    plot(r,out(1,:),':',r,out(2,:));
    title('estimated defocus vs actual');
    xlabel('actual');
    ylabel('estimated');
    figure;
    plot(r,out(3,:),r,out(4,:))
    title('minimum values for actual and estimate');
else
    plot(r,(out(1,:)),r,out(4,:),':');
    title('estimated defocus vs actual');
    xlabel('actual');
    ylabel('estimated');
    figure;
    plot(r,(out(2,:)),r,out(5,:),':');
    title('estimated aon vs actual');
    xlabel('actual');
    ylabel('estimated');
    figure;
    plot(r,(out(3,:)),r,out(6,:),':');
    title('estimated adiag vs actual');
    xlabel('actual');
    ylabel('estimated');
end
% set new WD/Stig from algorithm
% z=O(1);
% a_on=O(2);
% a_diag=O(3);