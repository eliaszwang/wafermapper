%MAPFoSt test using simulated images
close all
clear
raw=load('D:\Academics\Research\Seung Research\test images\[0 0 0].mat');
focused=raw.I1;
focused2=raw.I2;
% raw=(raw.out);
% raw1=(raw(:,1:1024));
% raw2=(raw(:,1025:2048));
% raw3=raw(1:512,1:512);


single=1;
tic;
out=[];
r=1:10; %loop over different (actual) aberrations
for i=r
raw=load(['D:\Academics\Research\Seung Research\test images\[' num2str(i) ' 7 -7].mat']);
%% Initialize/calculate constants
%Calculate fft of two images
I1=double(raw.I1);
I2=double(raw.I2);
fI1=fft2(I1); %image should have dimension 2^n for faster FFT
fI2=fft2(I2);
height=size(fI1,1);
width=size(fI1,2);
FOV=8.511;
Acc=5;
PixSize = FOV/height; % um per pixel
A_max=80; %set max defocus and astigmatism to 80um-based of paper, needs to be changed
A_sigma=80;
NA= 0.752 / (PixSize* (Acc*1000)^0.5);
sigma=5; %estimated Gaussian noise (approximation for shot noise), rad/um (maybe calculate later)
%[Kx, Ky]=meshgrid((mod(0.5+[0:width-1]/width,1)-0.5)*(6.28/FOV),(mod(0.5+[0:height-1]/height,1)-0.5)*(6.28/FOV)); %units are rad/um?
[Kx, Ky]=meshgrid((circshift([0:width-1]/width,width/2,2)-0.5)*(6.28/FOV),(circshift([0:width-1]/width,width/2,2)-0.5)*(6.28/FOV)); %units are rad/um?

if single
    %test initial aberration
    A=raw.A;
    % hardcoded test aberrations (defocus only)
    T1=raw.T1+i-1; %defocus in [um]
    T2=raw.T2+i-1;
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




% noise=fft2(normrnd(0,sigma,height,width));
% fI1=fI1./MTF(Kx,Ky,A+T1);
noise=fft2(normrnd(0,sigma,height,width));
fI3=fft2(focused).*MTF(Kx,Ky,10)+noise; %simulated image
I1=(ifft2(fI1));
I2=(ifft2(fI2));
I3=(ifft2(fI3));


% figure;
% imshow([fftshift(fI1) fftshift(fI2)]);
% figure;
% imshow([uint8(255*I1/max(max(I1))) uint8(255*I2/max(max(I2))) uint8(255*I3/max(max(I3)))]);
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
% MAPout=[];
% temp=-10:10;
% for j=temp
%     MAPout=[MAPout MAP(j,I1,I2,T1,T2,FOV,Acc,1)];
% end
% h=figure;
% plot(temp,MAPout);
% title(['simulation -ln(P(A)) vs A for A=' num2str(A)]);
%saveas(h,['D:\Academics\Research\Seung Research\Analysis plots\simulation raw3 -ln(P(A)) vs A for A=' num2str(A) '.jpg']);

end
toc;

if single
    figure;
    plot(r,out(1,:),':',r,out(2,:),'*');
    title('estimated defocus vs actual');
    xlabel('actual');
    ylabel('estimated');
    figure;
    plot(r,out(3,:),r,out(4,:))
    title('minimum values for actual and estimate');
else
    plot(r,(out(1,:)),':',r,out(4,:));
    title('estimated defocus vs actual');
    xlabel('actual');
    ylabel('estimated');
    figure;
    plot(r,(out(2,:)),':',r,out(5,:));
    title('estimated aon vs actual');
    xlabel('actual');
    ylabel('estimated');
    figure;
    plot(r,(out(3,:)),':',r,out(6,:));
    title('estimated adiag vs actual');
    xlabel('actual');
    ylabel('estimated');
end
% 
% %MAP function to maximize with respect to aberration vector, A
% %T1 and T2 are the known test aberrations
% gaussian =1;
% %% Initialize/calculate constants
% %Calculate fft of two images
% fI1=fft2(I1); %image should have dimension 2^n for faster FFT
% fI2=fft2(I2);
% height=size(I1,1);
% width=size(I1,2);
% PixSize = FOV/height; % um per pixel
% A_max=80; %set max defocus and astigmatism to 80um-based of paper, needs to be changed
% A_sigma=80; %sigma for gaussian prior
% %NA=5; %set NA to 5 mrad-based of paper, needs to be changed
% NA= 0.752 / (PixSize* (Acc*1000)^0.5); %confirm units for PixSize
% %sigma=10; %estimated Gaussian noise (approximation for shot noise), rad/um (maybe calculate later)
% sigma =mean([std(double(I1(:))), std(double(I2(:)))]);
% cutoffx=int32(floor(0.125*width)); %cutoff for k's used based on 25% k_nyquist, cycles/pixel
% cutoffy=int32(floor(0.125*height)); %use for selecting subset of I/K e.g. K([1:cutoffy end+1-cutoffy:end],[1:cutoffx end+1-cutoffx:end])
% [Kx, Ky]=meshgrid((mod(0.5+[0:width-1]/width,1)-0.5)*(6.28/FOV),(mod(0.5+[0:height-1]/height,1)-0.5)*(6.28/FOV));
% % determine which wave vectors to use
% % Kx=Kx([1:cutoffy end+1-cutoffy:end],[1:cutoffx end+1-cutoffx:end]);
% % Ky=Ky([1:cutoffy end+1-cutoffy:end],[1:cutoffx end+1-cutoffx:end]);
% % fI1=fI1([1:cutoffy end+1-cutoffy:end],[1:cutoffx end+1-cutoffx:end]);
% % fI2=fI2([1:cutoffy end+1-cutoffy:end],[1:cutoffx end+1-cutoffx:end]);
% Kx=Kx(1,:);
% Ky=Ky(1,:);
% fI1=fI1(1,:);
% fI2=fI2(1,:);
% %precompute some matrices
% Kx2=Kx.^2;
% Ky2=Ky.^2;
% 
% if single
%     %for single aberration mode
%     Kx2Ky2=Kx2+Ky2;
%     p_A=@(A)  (max(abs(A))<=A_max)/(2*A_max); %uniform distribtion over cube of side 2*A_max
%     MTF=@(Kx,Ky,A) exp(-0.125*(NA^2)*(Kx2Ky2)*A^2);
%     dMTFdz=@(Kx,Ky,A) -0.25*(NA^2)*(Kx2Ky2)*A.*exp(-0.125*(NA^2)*(Kx2Ky2)*A^2);
%     MTF1=MTF(Kx,Ky,A+T1);
%     MTF2=MTF(Kx,Ky,A+T2);
%     df=sum(sum( (2*(MTF1.*dMTFdz(Kx,Ky,A+T2)-MTF2.*dMTFdz(Kx,Ky,A+T1)).*(fI2.*MTF2+fI1.*MTF1).*(fI1.*MTF2-fI2.*MTF1))./(2*sigma^2*(MTF1.^2+MTF2.^2).^2) ));
%     if gaussian
%         p_A=@(A)  (1/(2.5066*A_sigma))*exp(-A^2/(2*A_sigma^2)); %gaussian prior
%         df=df+A/A_sigma^2;
%     end
%     f=-log(p_A(A))+sum(sum( ((fI2.*MTF1 - fI1.*MTF2).^2) ./ (2*sigma^2*(MTF1.^2+MTF2.^2)) ));
% else
%     %% Setup probability functions for MAPFoSt
%     p_A=@(A)  (max([abs(A(1)),abs(A(2)),abs(A(3))])<=A_max)/((2*A_max)^3); %uniform distribtion over cube of side 2*A_max
%     MTF=@(Kx,Ky,A) exp(-0.125*(NA^2)*(2*A(2)*(Kx2 - Ky2)*A(1) - A(3)*Kx.*Ky*A(1) + (Kx2+Ky2)*(A(2)^2 + A(3)^2 + A(1)^2)));
%     MTF1=MTF(Kx,Ky,A+T1);
%     MTF2=MTF(Kx,Ky,A+T2);
%     %MTF equation used to calculate partial derivatives in Wolfram
%     %exp(-0.125*(N^2)*(2a*(x^2-y^2)*z-b*x*y*z+(x^2+y^2)(a^2+b^2+z^2)))
%     dMTFd1=@(Kx,Ky,A) -0.125*(NA^2)*(2*A(2)*(Kx2 - Ky2) + 2*A(1)*(Kx2 + Ky2) - A(3)*Kx.*Ky).*exp(-0.125*(NA^2)*(2*A(2)*(Kx2 - Ky2)*A(1) - A(3)*Kx.*Ky*A(1) + (Kx2+Ky2)*(A(2)^2 + A(3)^2 + A(1)^2)));
%     dMTFd2=@(Kx,Ky,A) -0.125*(NA^2)*(2*A(2)*(Kx2 + Ky2) + 2*A(1)*(Kx2 - Ky2)).*exp(-0.125*(NA^2)*(2*A(2)*(Kx2 - Ky2)*A(1) - A(3)*Kx.*Ky*A(1) + (Kx2+Ky2)*(A(2)^2 + A(3)^2 + A(1)^2)));
%     dMTFd3=@(Kx,Ky,A) -0.125*(NA^2)*(2*A(3)*(Kx2 + Ky2) - A(1)*Kx.*Ky).*exp(-0.125*(NA^2)*(2*A(2)*(Kx.^2 - Ky.^2)*A(1) - A(3)*Kx.*Ky*A(1) + (Kx2+Ky2)*(A(2)^2 + A(3)^2 + A(1)^2)));
% 
%     % d/dx (a*f(x+i)-b*f(x+j))^2/(c*(f(x+i)^2+f(x+j)^2))
%     df(1)= sum(sum( (2*(MTF1.*dMTFd1(Kx,Ky,A+T2)-MTF2.*dMTFd1(Kx,Ky,A+T1)).*(fI2.*MTF2+fI1.*MTF1).*(fI1.*MTF2-fI2.*MTF1))./(2*sigma^2*(MTF1.^2+MTF2.^2).^2) ));
%     df(2)= sum(sum( (2*(MTF1.*dMTFd2(Kx,Ky,A+T2)-MTF2.*dMTFd2(Kx,Ky,A+T1)).*(fI2.*MTF2+fI1.*MTF1).*(fI1.*MTF2-fI2.*MTF1))./(2*sigma^2*(MTF1.^2+MTF2.^2).^2) ));
%     df(3)= sum(sum( (2*(MTF1.*dMTFd3(Kx,Ky,A+T2)-MTF2.*dMTFd3(Kx,Ky,A+T1)).*(fI2.*MTF2+fI1.*MTF1).*(fI1.*MTF2-fI2.*MTF1))./(2*sigma^2*(MTF1.^2+MTF2.^2).^2) ));
%     f=-log(p_A(A))+sum(sum( ((fI2.*MTF1 - fI1.*MTF2).^2) ./ (2*sigma^2*(MTF1.^2+MTF2.^2)) ));
%     %f=sum(sum( ((fI2.*MTF1 - fI1.*MTF2).^2) ./ (2*sigma^2*(MTF1.^2+MTF2.^2)) ));
% end


%%
clear;
v=VideoWriter('D:\Academics\Research\Seung Research\test.avi');
v.FrameRate=1;
open(v);
x=0:40;
for i=x
raw=load(['D:\Academics\Research\Seung Research\MAPFoSt-test-images\test images 6_28_16\[' num2str(i) ' 7 -7].mat']);
FOV=8.511;
I1=double(raw.I1);
I2=double(raw.I2);
fI1=fftshift(fft2(I1)); %image should have dimension 2^n for faster FFT
fI2=fftshift(fft2(I2));
height=size(fI1,1);
width=size(fI1,2);
[Kx, Ky]=meshgrid(([0:width-1]/width-0.5)*(6.28/FOV),([0:height-1]/height-0.5)*(6.28/FOV)); %units are rad/um?
r=abs(-0.125*1.2795^2*(Kx.^2+Ky.^2).*(2*raw.A*(raw.T1-raw.T2)+raw.T1^2-raw.T2^2));
l=abs(log(fI1./fI2));
%imshow(l/max(l(:)));
imshow(l);
writeVideo(v,frame2im(getframe(gcf)));
end
close(v);