function [MSE]=MAPFoSt_NA_helper(NA,imagepath,imagedir)
    r=1:length(imagedir);
    out=[];
    for i=r
        raw=load([imagepath imagedir(i).name]);
        I1=double(raw.I1);
        I2=double(raw.I2);
        %subsample image
        I1=I1(1:4:1024,1:4:1024);
        I2=I2(1:4:1024,1:4:1024);

        %% Initialize/calculate constants
        %Calculate fft of two images
        fI1=fft2(I1); %image should have dimension 2^n for faster FFT
        fI2=fft2(I2);
        height=size(I1,1);
        width=size(I1,2);
        FOV=8.511;
        sigma =mean([std(double(I1(:))), std(double(I2(:)))]); %sigma for real space, approximation for shot noise
        %[Kx, Ky]=meshgrid((mod(0.5+[0:width-1]/width,1)-0.5)*(6.28/FOV),(mod(0.5+[0:height-1]/height,1)-0.5)*(6.28/FOV)); % use mod instead of cirshift for backwards compatibilty, units are rad/um
        [Kx, Ky]=meshgrid((circshift([0:width-1]/width,width/2,2)-0.5)*(6.28*width/FOV),(circshift([0:height-1]/height,height/2,2)-0.5)*(6.28*height/FOV)); %units are rad/um

        %test initial aberration
        A=raw.A(1);
        % hardcoded test aberrations (defocus only)
        T1=raw.T1(1); %defocus in [um]
        T2=raw.T2(1);
        init=2;



        p.length=20;
        p.method='BFGS';
        p.verbosity=0;
        p.MFEPLS = 30;   % Max Func Evals Per Line Search
        p.MSR = 100;                % Max Slope Ratio default
        O=real(minimize(init,@MAP,p,fI1,fI2,T1,T2,NA,sigma,Kx,Ky));
        out=[out [A';O';MAP(A,fI1,fI2,T1,T2,NA,sigma,Kx,Ky);MAP(O,fI1,fI2,T1,T2,NA,sigma,Kx,Ky)]];

    end
    MSE=[immse(out(1,:),out(2,:))];
end