function [f, df]=MAP(A,fI1,fI2,T1,T2,NA,sigma,Kx, Ky)
%   MAP function to maximize with respect to aberration vector, A
%   fI1,fI2: fft2 of the images
%   T1, T2: known test aberrations
%   NA: numerical aperature (empirically estimated)
%   sigma: image noise std
%   Kx,Ky: wave vector (meshgrid)

gaussian =0; %gaussian prior flag

%% Initialize/calculate constants
height=size(fI1,1);
width=size(fI1,2);
A_max=50; %set max defocus to 50um
A_sigma=20; %sigma for gaussian prior
NA2=NA^2;
cutoffx=int32(floor(0.25*width)); %cutoff for k's used based on 50% k_nyquist, cycles/pixel
cutoffy=int32(floor(0.25*height)); %use for selecting subset of I/K e.g. K([1:cutoffy end+1-cutoffy:end],[1:cutoffx end+1-cutoffx:end])
%% select which wave vectors to use (based on cutoff)
Kx=Kx([1:cutoffy end+1-cutoffy:end],[1:cutoffx end+1-cutoffx:end]);
Ky=Ky([1:cutoffy end+1-cutoffy:end],[1:cutoffx end+1-cutoffx:end]);
fI1=fI1([1:cutoffy end+1-cutoffy:end],[1:cutoffx end+1-cutoffx:end]);
fI2=fI2([1:cutoffy end+1-cutoffy:end],[1:cutoffx end+1-cutoffx:end]);
%% precompute some matrices
Kx2=Kx.^2;
Ky2=Ky.^2;
Kx2Ky2=Kx2+Ky2;
p_A=@(A)  (max(abs(A))<=A_max)/(2*A_max); % define uniform distribtion over line of length 2*A_max
MTF=@(Kx,Ky,A) exp(-0.125*(NA2)*(Kx2Ky2)*A^2); % define gaussian MTF
dMTFdz=@(Kx,Ky,A) -0.25*(NA2)*(Kx2Ky2)*A.*exp(-0.125*(NA2)*(Kx2Ky2)*A^2); % define derivative of MTF
MTF1=MTF(Kx,Ky,A+T1);
MTF2=MTF(Kx,Ky,A+T2);
MTF12=MTF1.^2;
MTF22=MTF2.^2;

likelihood=((fI2.*MTF1 - fI1.*MTF2).^2) ./ (2*sigma^2*(MTF12+MTF22)+1e-20); %this is the matrix of log likelihoods (ie the exponentiated part of eq 27)

df=(2*(MTF1.*dMTFdz(Kx,Ky,A+T2)-MTF2.*dMTFdz(Kx,Ky,A+T1)).*(fI2.*MTF2+fI1.*MTF1).*(fI1.*MTF2-fI2.*MTF1))./(2*sigma^2*(MTF12+MTF22).^2+1e-20) ;
df=sum(sum( df(likelihood>0) )); %only take the df terms from the ones we're using to calculate f
if gaussian
    p_A=@(A)  (1/(2.5066*A_sigma))*exp(-A^2/(2*A_sigma^2)); % define gaussian prior
    df=df+A/A_sigma^2;
end
f=-log(p_A(A))+sum(sum( likelihood(likelihood>0) )); %take positive log likelihoods
end