function [Imf] = mexHatSection2(I,s1,s2)

%% process I
kSize = [256 256];


kRad = (kSize + 1)/2;
kern = zeros(kSize);

[y x z] = ind2sub(kSize,find(kern==0));
dists = sqrt(((y-kRad(1))).^2 + ((x - kRad(2))).^2);

cKern = 1 * exp(-.5 * (dists/s1).^2);
cKern = cKern/sum(cKern(:));
sKern = 1 * exp(-.5 * (dists/s2).^2);
sKern = sKern/sum(sKern(:));
kern(:) = cKern - sKern;

% subplot(2,1,1)
% figure(3);
% clf;
% plot(kern(round(kRad(1)),:))

%% Convolve

Itemp = fastCon(I,kern);
pixClip = s1;
Imf  = Itemp * 0;
Imf(pixClip +1:end-pixClip,pixClip +1:end-pixClip)=Itemp(pixClip +1:end-pixClip,pixClip +1:end-pixClip);


themean=mean(Imf(:));
thestd=std(Imf(:));
Imf=(Imf-themean)*(128.0/(3.5*thestd));
Imf(Imf<-128)=-128;
Imf(Imf>127)=127;

Imf=Imf+128;

%Imf(Imf<0)=0;


%Imf = Imf*(256/max(Imf(:)));
Imf = 255-Imf;
