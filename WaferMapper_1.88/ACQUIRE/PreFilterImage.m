function [I] = PreFilterImage(I,Options)
if ~exist('Options','var')
    Options.MexicanHat
    Options.Inner=1;
    Options.Outer=50;
    Options.ReNormalize=1; 
    Options.Sigma=3.5;
end
if ~isfield(Options,'MexicanHat')
    Options.MexicanHat=0;
end
if ~isfield(Options,'Inner')
    Options.Inner=50;
end
if ~isfield(Options,'Outer')
    Options.Outer=1;
end
if ~isfield(Options,'ReNormalize')
    Options.ReNormalize=0;
end
if ~isfield(Options,'Sigma')
    Options.Sigma=3.5;
end
%% process I
if Options.MexicanHat
    kSize = [3*Options.Outer 3*Options.Outer];


    kRad = (kSize + 1)/2;
    kern = zeros(kSize);

    [y x z] = ind2sub(kSize,find(kern==0));
    dists = sqrt(((y-kRad(1))).^2 + ((x - kRad(2))).^2);

    cKern = 1 * exp(-.5 * (dists/Options.Outer).^2);
    cKern = cKern/sum(cKern(:));
    sKern = 1 * exp(-.5 * (dists/Options.Inner).^2);
    sKern = sKern/sum(sKern(:));
    kern(:) = cKern - sKern;

    % subplot(2,1,1)
    % figure(3);
    % clf;
    % plot(kern(round(kRad(1)),:))

    %% Convolve

    Itemp = fastCon(I,kern);
    pixClip = Options.Outer;
    Imf  = Itemp * 0;
    Imf(pixClip +1:end-pixClip,pixClip +1:end-pixClip)=Itemp(pixClip +1:end-pixClip,pixClip +1:end-pixClip);
    I=Imf;
  
end
if Options.ReNormalize
    themean=mean(I(:));
    thestd=std(I(:));
    I=(I-themean)*(128.0/(Options.Sigma*thestd));
    I(I<-128)=-128;
    I(I>127)=127;
    I=I+128;
end
if Options.MexicanHat
      I = 255-I;
end
%Imf(Imf<0)=0;


%Imf = Imf*(256/max(Imf(:)));

