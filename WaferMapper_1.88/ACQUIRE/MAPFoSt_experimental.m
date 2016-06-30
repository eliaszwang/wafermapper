%MAPFoSt test using real images
close all
clear
raw=load('D:\Academics\Research\Seung Research\MAPFoSt-test-images\test images 6_22_16\[0 0 0]');
focused=raw.I1; %focused image for reference
focused2=raw.I2;
single=1;
tic;
out=[];
r=0:80;
for i=r
raw=load(['D:\Academics\Research\Seung Research\MAPFoSt-test-images\test images 6_28_16\[' num2str(i) ' 15 -15].mat']);
I1=double(raw.I1);
I2=double(raw.I2);
% figure;
% subplot(2,1,1);
% imhist(I1);
% subplot(2,1,2);
% imhist(I2);

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

if single
    %test initial aberration
    A=raw.A;
    % hardcoded test aberrations (defocus only)
    T1=raw.T1; %defocus in [um]
    T2=raw.T2;
    init=2;
else
    %test initial aberration
    A=[raw.A 0 0];
    % hardcoded test aberrations (defocus only)
    T1=[raw.T1 0 0]; %defocus in [um]
    T2=[raw.T2 0 0];
    init=[0 0 0];
end

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

% fOmap=(fI1.*MTF(Kx,Ky,A+T1)+fI2.*MTF(Kx,Ky,A+T2))./(MTF(Kx,Ky,A+T1).^2+MTF(Kx,Ky,A+T2).^2);
% Omap=ifft2(fOmap);
% Omap=uint8(255*Omap/max(max(Omap)));
% figure;
% imshow(Omap);
% disp(['MSE ' num2str(A) ': ' num2str(immse(Omap,focused))]);

p.length=20;
p.method='BFGS';
p.verbosity=0;
p.MFEPLS = 10;   % Max Func Evals Per Line Search
p.MSR = 100;                % Max Slope Ratio default
O=minimize(init,@MAP,p,fI1,fI2,T1,T2,NA,sigma,Kx,Ky,single);
out=[out [A';O';MAP(A,fI1,fI2,T1,T2,NA,sigma,Kx,Ky,single);MAP(O,fI1,fI2,T1,T2,NA,sigma,Kx,Ky,single)]];
% if single
%     % plot 1-D MAP
%     fout=[];
%     dfout=[];
%     temp=-10:20;
%     for j=temp
%         [f,df]=MAP(j,I1,I2,T1,T2,FOV,Acc,single);
%         fout=[fout f ];
%         dfout=[dfout df];
%     end
%     h=figure;
%     plot(temp,fout);
%     yyaxis right;
%     plot(temp,dfout);
%     title(['experimental -ln(P(A)) vs A for A=' num2str(A)]);
%     saveas(h,['D:\Academics\Research\Seung Research\Analysis plots\experimental second set 15 -ln(P(A)) vs A for A=' num2str(A) '.jpg']);
% end
end
toc;

if single
    ind=(out(3,:)-out(4,:)>=0);
    ind2=(out(3,:)-out(4,:)<0);
    figure;
    plot(r(ind),out(1,ind),':',r(ind),out(2,ind),'*',r(ind2),out(2,ind2),'+');
    title('estimated defocus vs actual');
    xlabel('actual');
    ylabel('estimated');
    figure;
    plot(r,out(3,:),r,out(4,:),'*')
    title('minimum values for actual and estimate');
else
    plot(r,(out(1,:)),':',r,out(4,:),'*');
    title('estimated defocus vs actual');
    xlabel('actual');
    ylabel('estimated');
    figure;
    plot(r,(out(2,:)),':',r,out(5,:),'*');
    title('estimated aon vs actual');
    xlabel('actual');
    ylabel('estimated');
    figure;
    plot(r,(out(3,:)),':',r,out(6,:),'*');
    title('estimated adiag vs actual');
    xlabel('actual');
    ylabel('estimated');
end
% set new WD/Stig from algorithm
% z=O(1);
% a_on=O(2);
% a_diag=O(3);