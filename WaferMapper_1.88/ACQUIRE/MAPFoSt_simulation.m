%MAPFoSt test using simulated images
close all
clear
raw=load('../../../MAPFoSt-test-images/test images 6_22_16/[0 0 0].mat');
% raw=(raw.out);
% raw1=(raw(:,1:1024));
% raw2=(raw(:,1025:2048));
% raw3=raw(1:512,1:512);
% psf=zeros(1024,1024);
% psf(500:525,500:525)=1;
% psf=psf/sum(sum(psf));
% mtf=fft2(psf);


single=0;
tic;
out=[];
r=0:20; %loop over different (actual) aberrations
for i=r
%% Initialize/calculate constants
I1=raw.I1;
I2=raw.I1;
% I1=filter2(fspecial('average',3),raw.I1)/255;
% I2=filter2(fspecial('average',3),raw.I2)/255;
%subsample image
I1=I1(1:4:1024,1:4:1024);
I2=I2(1:4:1024,1:4:1024);
% t=256;
% I1=I1(1:t,1:t);
% I2=I2(1:t,1:t);
%Calculate fft of two images
fI1=fft2(double(I2)); %image should have dimension 2^n for faster FFT
fI2=fft2(double(I2));
height=size(I1,1);
width=size(I1,2);
FOV=8.511;
Acc=5;
PixSize = FOV/height; % um per pixel
A_max=80; %set max defocus and astigmatism to 80um-based of paper, needs to be changed
%NA= sqrt(40)*0.752 / (PixSize* (Acc*1000)^0.5);
NA= 0.5596*height / ((Acc*1000)^0.5); %empirically determined constant;
sigmaI=1; %estimated Gaussian noise (approximation for shot noise), rad/um (maybe calculate later)
sigma =mean([std(double(I1(:))), std(double(I2(:)))]); %sigma for real space
%[Kx, Ky]=meshgrid((mod(0.5+[0:width-1]/width,1)-0.5)*(6.28/FOV),(mod(0.5+[0:height-1]/height,1)-0.5)*(6.28/FOV)); %units are rad/um?
[Kx, Ky]=meshgrid((circshift([0:width-1]/width,width/2,2)-0.5)*(6.28/FOV),(circshift([0:height-1]/height,height/2,2)-0.5)*(6.28/FOV)); %units are rad/um?

if single
    %test initial aberration
    A=i;
    % hardcoded test aberrations (defocus only)
    T1=15; %defocus in [um]
    T2=-15;
    init=2;
    MTF=@(Kx,Ky,A) exp(-0.125*(NA^2)*(Kx.^2+Ky.^2)*A^2);
else
    %test initial aberration
    A=[i -2 3];
    % hardcoded test aberrations (defocus only)
    T1=[15 0 0]; %defocus in [um]
    T2=[-15 0 0];
    init=[2 2 2];
    MTF=@(Kx,Ky,A) exp(-0.125*(NA^2)*(-2*A(2)*(Kx.^2 - Ky.^2)*A(1) - 4*A(3)*Kx.*Ky*A(1) + (Kx.^2+Ky.^2)*(A(2)^2 + A(3)^2 + A(1)^2)));
end




noise=fft2(normrnd(0,sigmaI,height,width));
fI1=fI1.*MTF(Kx,Ky,A+T1)+noise;
noise=fft2(normrnd(0,sigmaI,height,width));
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
% checkgrad('MAP', randn(1,1), 1e-5,fI1,fI2,T1',T2',NA,sigma,Kx,Ky,single);
p.length=20;
p.method='BFGS';
p.verbosity=1;
p.MFEPLS = 30;   % Max Func Evals Per Line Search
p.MSR = 100;                % Max Slope Ratio default
O=real(minimize(init,@MAP,p,fI1,fI2,T1,T2,NA,sigma,Kx,Ky,single));
if ~single
    O=real(minimize([(O(1)) 0 0],@MAP,p,fI1,fI2,T1,T2,NA,sigma,Kx,Ky,single));
    %O=real(minimize([O(1) O(2) 0],@MAP,p,fI1,fI2,T1,T2,NA,sigma,Kx,Ky,single));
end
out=[out [A';O';MAP(A,fI1,fI2,T1,T2,NA,sigma,Kx,Ky,single);MAP(O,fI1,fI2,T1,T2,NA,sigma,Kx,Ky,single)]];
% if single
%     %plot 1-D MAP
%     fout=[];
%     dfout=[];
%     temp=-10:10;
%     for j=temp
%         [f,df]=MAP(j,I1,I2,T1,T2,FOV,Acc,single);
%         fout=[fout f ];
%         dfout=[dfout df];
%     end
%     h=figure;
%     plot(temp,fout);
%     yyaxis right;
%     plot(temp,dfout);
%     title(['simulation -ln(P(A)) vs A for A=' num2str(A)]);
%     %saveas(h,['D:\Academics\Research\Seung Research\Analysis plots\simulation no abs -ln(P(A)) vs A for A=' num2str(A) '.jpg']);
% end


end
% MSE=immse(out(1,:),out(2,:));
% disp(['MSE: ' num2str(MSE)]);
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
    %R=[-0.1140 0.0111;0.0602 0.0332];
    R=[-0.1245 0.0093;0.0584 0.0323];
    out(9:10,:)=R*out(2:3,:);
    out(11:12,:)=R*out(5:6,:);
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
    plot(r,(out(9,:)),':',r(ind),out(11,ind),'*',r(ind2),out(11,ind2),'+');
    title('estimated ax vs actual');
    xlabel('actual');
    ylabel('estimated');
    figure;
    plot(r,(out(10,:)),':',r(ind),out(12,ind),'*',r(ind2),out(12,ind2),'+');
    title('estimated ay vs actual');
    xlabel('actual');
    ylabel('estimated');
end
% set new WD/Stig from algorithm
% z=O(1);
% a_on=O(2);
% a_diag=O(3);