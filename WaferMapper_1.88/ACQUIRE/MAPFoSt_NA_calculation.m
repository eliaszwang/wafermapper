%Script for determining MAPFoSt NA based on experimental images
close all
clear
imagepath='../../../MAPFoSt-test-images/test images 6_28_16/';
imagedir=dir([imagepath '[* 15 -15].mat']);
imagedir=imagedir(1:20,1);
%% initial log scale search over 4 orders of magnitude
min_value=0.01; %initialization, NA must be within 2 orders of magnitude of this value
tic;
NA_range=logspace(log10(min_value)-2,log10(min_value)+2,10); %NA values
MSE=[];
for NA=NA_range
    MSE=[MSE MAPFoSt_NA_helper(NA,imagepath,imagedir);];
end
toc;
figure;
plot(NA_range,MSE,'*');
xlabel('NA');
ylabel('MSE');
min_index=find(MSE==min(MSE));
min_value=NA_range(min_index);

%% log search within 2 orders of magnitude
tic;
NA_range=logspace(log10(min_value)-1,log10(min_value)+1,11); %NA values
MSE=[];
for NA=NA_range
    MSE=[MSE MAPFoSt_NA_helper(NA,imagepath,imagedir);];
end
toc;
figure;
plot(NA_range,MSE,'*');
xlabel('NA');
ylabel('MSE');
min_index=find(MSE==min(MSE));
min_value=NA_range(min_index);

%% linear search 
tic;
NA_range=linspace(min_value*10^-0.2,min_value*10^0.2,10); %NA values
MSE=[];
for NA=NA_range
    MSE=[MSE MAPFoSt_NA_helper(NA,imagepath,imagedir);];
end
toc;
figure;
plot(NA_range,MSE,'*');
xlabel('NA');
ylabel('MSE');
min_index=find(MSE==min(MSE));
min_value=NA_range(min_index);

%% finer linear search, if necessary

tic;
NA_range=linspace(0.9*min_value,1.1*min_value,20); %NA values
MSE=[];
for NA=NA_range
    MSE=[MSE MAPFoSt_NA_helper(NA,imagepath,imagedir);];
end
toc;
figure;
plot(NA_range,MSE,'*');
xlabel('NA');
ylabel('MSE');
min_index=find(MSE==min(MSE));
min_value=NA_range(min_index);

%{
raw=load(['../../../MAPFoSt-test-images/test images 6_28_16/[' num2str(3) ' 15 -15].mat']);
out=[];
tic;
r=0.6:0.1:8;
for i=r
    I1=double(raw.I1);
I2=double(raw.I2);
%subsample image
I1=I1(1:4:1024,1:4:1024);
I2=I2(1:4:1024,1:4:1024);
% I1=I1(257:768,257:768);
% I2=I2(257:768,257:768);
% I1=I1(385:640,385:640);
% I2=I2(385:640,385:640);
% t=256;
% I1=I1(1:t,1:t);
% I2=I2(1:t,1:t);
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
NA=0.0079*i; %empirically determined NA
sigma =mean([std(double(I1(:))), std(double(I2(:)))]); %sigma for real space, approximation for shot noise
%[Kx, Ky]=meshgrid((mod(0.5+[0:width-1]/width,1)-0.5)*(6.28/FOV),(mod(0.5+[0:height-1]/height,1)-0.5)*(6.28/FOV)); % use mod instead of cirshift for backwards compatibilty, units are rad/um
[Kx, Ky]=meshgrid((circshift([0:width-1]/width,width/2,2)-0.5)*(6.28*width/FOV),(circshift([0:height-1]/height,height/2,2)-0.5)*(6.28*height/FOV)); %units are rad/um


    %test initial aberration
    A=raw.A(1);
    % hardcoded test aberrations (defocus only)
    T1=raw.T1(1); %defocus in [um]
    T2=raw.T2(1);
    init=2;

% close all
% figure;
% imshow([fftshift(fI1) fftshift(fI2)]);
% figure;
% imshow([uint8(255*I1/max(max(I1))) uint8(255*I2/max(max(I2)))]);
% 
% fOmap=(fI1.*MTF(Kx,Ky,A+T1)+fI2.*MTF(Kx,Ky,A+T2))./(MTF(Kx,Ky,A+T1).^2+MTF(Kx,Ky,A+T2).^2);
% Omap=ifft2(fOmap);
% Omap=(Omap/max(Omap(:)));
% figure;
% imshow(Omap);
% 
% MAP(A,fI1,fI2,T1,T2,NA,sigma,Kx, Ky,single);
% checkgrad('MAP', randn(1,1), 1e-5,fI1,fI2,T1',T2',NA,sigma,Kx,Ky,single);

% fOmap=(fI1.*MTF(Kx,Ky,A+T1)+fI2.*MTF(Kx,Ky,A+T2))./(MTF(Kx,Ky,A+T1).^2+MTF(Kx,Ky,A+T2).^2);
% Omap=ifft2(fOmap);
% Omap=uint8(255*Omap/max(max(Omap)));
% figure;
% imshow(Omap);
% disp(['MSE ' num2str(A) ': ' num2str(immse(Omap,focused))]);

p.length=20;
p.method='BFGS';
p.verbosity=0;
p.MFEPLS = 30;   % Max Func Evals Per Line Search
p.MSR = 100;                % Max Slope Ratio default
O=real(minimize(init,@MAP,p,fI1,fI2,T1,T2,NA,sigma,Kx,Ky));

out=[out [A';O';MAP(A,fI1,fI2,T1,T2,NA,sigma,Kx,Ky);MAP(O,fI1,fI2,T1,T2,NA,sigma,Kx,Ky)]];

% if single
%     % plot 1-D MAP
%     fout=[];
%     dfout=[];
%     temp=-10:20;
%     for j=temp
%         [f,df]=MAP(j,fI1,fI2,T1,T2,NA,sigma,Kx,Ky,single);
%         fout=[fout f ];
%         dfout=[dfout df];
%     end
%     h=figure;
%     plot(temp,fout);
%     yyaxis right;
%     plot(temp,dfout);
%     title(['experimental -ln(P(A)) vs A for A=' num2str(A)]);
%     %saveas(h,['D:\Academics\Research\Seung Research\Analysis plots\experimental second set 15 -ln(P(A)) vs A for A=' num2str(A) '.jpg']);
% else
%     %plot 2-D aon-adiag
%     fout=zeros(41,41);
%     for j=-10:0.5:10
%         for k=-10:0.5:10
%             fout(uint8(j*2+21),uint8(k*2+21))=MAP([A(1) j k],fI1,fI2,T1,T2,NA,sigma,Kx,Ky,single);
%         end
%     end
%     figure;
%     imagesc(real(fout));
% end
end
toc;
close all
figure;
plot(r,out(2,:));
xlabel('NA')
ylabel('estimation')
%}