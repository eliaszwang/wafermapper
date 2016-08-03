%MAPFoSt test using real images
close all
%clear
raw=load('../../../MAPFoSt-test-images/test images 6_22_16/[0 0 0].mat');
focused=raw.I1; %focused image for reference
focused2=raw.I2;
single=0;
tic;
out=[];
r=-10:10;
for i=r
%raw=load(['../../../MAPFoSt-test-images/test images 7_21_16/defocusy[' num2str(i) ' 0 1][15 0 0][-15 0 0]PixSize8.mat']);
raw=load(['../../../MAPFoSt-test-images/test images 7_21_16/stigx[0 ' num2str(i) ' 0][15 0 0][-15 0 0]PixSize8.mat']);
%raw=load(['../../../MAPFoSt-test-images/test images 6_28_16/[' num2str(i) ' 15 -15].mat']);
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
Acc=5;
PixSize = FOV/height; % um per pixel
% NA= sqrt(40)*0.752 / (8.511/height* (Acc*1000)^0.5);
NA= 0.5596*height / ((Acc*1000)^0.5);
%sigma=5; %estimated Gaussian noise (approximation for shot noise), rad/um (maybe calculate later)
sigma =mean([std(double(I1(:))), std(double(I2(:)))]); %sigma for real space
%[Kx, Ky]=meshgrid((mod(0.5+[0:width-1]/width,1)-0.5)*(6.28/FOV),(mod(0.5+[0:height-1]/height,1)-0.5)*(6.28/FOV)); %units are rad/um?
[Kx, Ky]=meshgrid((circshift([0:width-1]/width,width/2,2)-0.5)*(6.28/FOV),(circshift([0:width-1]/width,width/2,2)-0.5)*(6.28/FOV)); %units are rad/um?

if single
    %test initial aberration
    A=raw.A(1);
    % hardcoded test aberrations (defocus only)
    T1=raw.T1(1); %defocus in [um]
    T2=raw.T2(1);
    init=2;
else
    %test initial aberration
    A=[raw.A ];
    % hardcoded test aberrations (defocus only)
    T1=[raw.T1 ]; %defocus in [um]
    T2=[raw.T2 ];
    init=[2 2 2];
end

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
O=real(minimize(init,@MAP,p,fI1,fI2,T1,T2,NA,sigma,Kx,Ky,single));
if ~single
    O=real(minimize([(O(1)) 2 2],@MAP,p,fI1,fI2,T1,T2,NA,sigma,Kx,Ky,single));
    %O=real(minimize([O(1) O(2) 0],@MAP,p,fI1,fI2,T1,T2,NA,sigma,Kx,Ky,single));
end
out=[out [A';O';MAP(A,fI1,fI2,T1,T2,NA,sigma,Kx,Ky,single);MAP(O,fI1,fI2,T1,T2,NA,sigma,Kx,Ky,single)]];

if single
    % plot 1-D MAP
    fout=[];
    dfout=[];
    temp=-10:20;
    for j=temp
        [f,df]=MAP(j,fI1,fI2,T1,T2,NA,sigma,Kx,Ky,single);
        fout=[fout f ];
        dfout=[dfout df];
    end
    h=figure;
    plot(temp,fout);
    yyaxis right;
    plot(temp,dfout);
    title(['experimental -ln(P(A)) vs A for A=' num2str(A)]);
    %saveas(h,['D:\Academics\Research\Seung Research\Analysis plots\experimental second set 15 -ln(P(A)) vs A for A=' num2str(A) '.jpg']);
else
    %plot 2-D aon-adiag
    fout=zeros(41,41);
    for j=-10:0.5:10
        for k=-10:0.5:10
            fout(uint8(j*2+21),uint8(k*2+21))=MAP([A(1) j k],fI1,fI2,T1,T2,NA,sigma,Kx,Ky,single);
        end
    end
    figure;
    imagesc(real(fout));
end
end
toc;
if single
    ind=(out(3,:)-out(4,:)>=0);
    ind2=(out(3,:)-out(4,:)<0);
    disp(num2str(immse(out(1,ind),out(2,ind))));
    figure;
    plot(r,out(1,:),':',r(ind),out(2,ind),'*',r(ind2),out(2,ind2),'+');
    title('estimated defocus vs actual');
    xlabel('actual');
    ylabel('estimated');
    figure;
    plot(r,out(3,:),r,out(4,:),'*')
    title('minimum values for actual and estimate');
else
    % rotation from aon,adiag to x,y
    %R=[-0.1245 0.0093;0.0584 0.0323]; %old MTF
    R=[0.1273 0.0436;-0.0548 0.1306]; %new MTF
    out(9:10,:)=R*out(5:6,:);
    ind=(out(7,:)-out(8,:)>=0);
    ind2=(out(7,:)-out(8,:)<0);
    figure;
    plot(r,(out(1,:)),':',r(ind),out(4,ind),'*',r(ind2),out(4,ind2),'+');
    title('estimated defocus vs actual');
    xlabel('actual');
    ylabel('estimated');
    figure;
    plot(r,(out(2,:)),':',r(ind),out(5,ind),'*',r(ind2),out(5,ind2),'+');
    title('estimated aon vs actual');
    xlabel('actual');
    ylabel('estimated');
    figure;
    plot(r,(out(3,:)),':',r(ind),out(6,ind),'*',r(ind2),out(6,ind2),'+');
    title('estimated adiag vs actual');
    xlabel('actual');
    ylabel('estimated');
    figure;
    plot(r,out(7,:),':',r,out(8,:),'*')
    title('minimum values for actual and estimate');
    figure;
    plot(r,(out(2,:)),':',r(ind),out(9,ind),'*',r(ind2),out(9,ind2),'+');
    title('estimated ax vs actual');
    xlabel('actual');
    ylabel('estimated');
    figure;
    plot(r,(out(3,:)),':',r(ind),out(10,ind),'*',r(ind2),out(10,ind2),'+');
    title('estimated ay vs actual');
    xlabel('actual');
    ylabel('estimated');
end
