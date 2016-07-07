function [z]=MAPFoSt(ImageHeightInPixels,ImageWidthInPixels,DwellTimeInMicroseconds,FileName,FOV)
% Elias Wang's implementation of MAPFoSt
% algorithm adapted from Jonas Binding (Low-Dosage Maximum-A-Posteriori Focusing and Stigmation)
% should replace "sm.Execute('CMD_AUTO_FOCUS_FINE')" calls
global GuiGlobalsStruct;
sm = GuiGlobalsStruct.MyCZEMAPIClass; %To shorten calls to global API variables in this function


% hardcoded aberrations
T1=15;
T2=-15;


frametime=sm.Get_ReturnTypeSingle('AP_FRAME_TIME')/1000;
frametime=0.3;

%begin MAPFoSt
CurrentWorkingDistance = sm.Get_ReturnTypeSingle('AP_WD');
disp('Beginning MAPFoSt...');
disp(['starting WD: ' num2str(10^6*CurrentWorkingDistance) 'um']);


% %*** START: This sequence is designed to release the SEM from Fibics control
%     sm.Execute('CMD_AUTO_FOCUS_FINE');
%     pause(0.5);
%     sm.Execute('CMD_ABORT_AUTO');
%     while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
%         pause(0.02);
%     end
%     sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance);
%     
%     pause(0.1);
   

%%Take first image
%implement errorcheck for test aberration
sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T1);
sm.Set_PassedTypeSingle('AP_SCANROTATION',0);
pause(frametime);
T1WD=sm.Get_ReturnTypeSingle('AP_WD');
disp(['first image WD: ' num2str(10^6*T1WD) 'um']);
sm.Fibics_WriteFOV(FOV);
%Wait for image to be acquired
sm.Fibics_AcquireImage(ImageWidthInPixels,ImageHeightInPixels,DwellTimeInMicroseconds,FileName);
while(sm.Fibics_IsBusy)
    pause(.01); %1
end

%Wait for file to be written
IsReadOK = false;
while ~IsReadOK
    IsReadOK = true;
    try
        I1 = imread(FileName);
    catch MyException
        IsReadOK = false;
        pause(0.1);
    end
end

delete(FileName);

% %*** START: This sequence is designed to release the SEM from Fibics control
%     sm.Execute('CMD_AUTO_FOCUS_FINE');
%     pause(0.5);
%     sm.Execute('CMD_ABORT_AUTO');
%     while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
%         pause(0.02);
%     end
%     sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance);
%     
%     pause(0.1);
%     


%%Take second image
%implement errorcheck for test aberration
sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T2);
sm.Set_PassedTypeSingle('AP_SCANROTATION',0);
pause(frametime);
T2WD=sm.Get_ReturnTypeSingle('AP_WD');
disp(['second image WD: ' num2str(10^6*T2WD) 'um']);
sm.Fibics_WriteFOV(FOV);
%Wait for image to be acquired
sm.Fibics_AcquireImage(ImageWidthInPixels,ImageHeightInPixels,DwellTimeInMicroseconds,FileName);
while(sm.Fibics_IsBusy)
    pause(.01); %1
end

%Wait for file to be written
IsReadOK = false;
while ~IsReadOK
    IsReadOK = true;
    try
        I2 = imread(FileName);
    catch MyException
        IsReadOK = false;
        pause(0.1);
    end
end

delete(FileName);

%*** START: This sequence is designed to release the SEM from Fibics control
    sm.Execute('CMD_AUTO_FOCUS_FINE');
    pause(0.5);
    sm.Execute('CMD_ABORT_AUTO');
    while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
        pause(0.02);
    end
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance);
    
    pause(0.1);

%get actual aberrations
    T1=10^6*(T1WD-CurrentWorkingDistance);
    T2=10^6*(T2WD-CurrentWorkingDistance);
    
%% MAPFoSt
%PixSize = FOV/ImageHeightInPixels; % um per pixel
Acc=5;
%NA= sqrt(40)*0.752 / (PixSize* (Acc*1000)^0.5);
NA= 0.5596*ImageHeightInPixels / ((Acc*1000)^0.5); %empirically determined constant;
sigma =mean([std(double(I1(:))), std(double(I2(:)))]); %sigma for real space
[Kx, Ky]=meshgrid((mod(0.5+[0:ImageWidthInPixels-1]/ImageWidthInPixels,1)-0.5)*(6.28/FOV),(mod(0.5+[0:ImageHeightInPixels-1]/ImageHeightInPixels,1)-0.5)*(6.28/FOV));
fI1=fft2(double(I1)); %image should have dimension 2^n for faster FFT
fI2=fft2(double(I2));
init=2;



p.length=20;
p.method='BFGS';
p.verbosity=1;
p.MFEPLS = 30;   % Max Func Evals Per Line Search
p.MSR = 100;     % Max Slope Ratio default
O=minimize(init,@MAP,p,fI1,fI2,T1,T2,NA,sigma,Kx,Ky,1);

%%
if max(abs(O))<80
    % set new WD/Stig from algorithm
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance-10^-6*real(O)); %change WD to testing defocus
    disp(['final WD: ' num2str(10^6*sm.Get_ReturnTypeSingle('AP_WD')) 'um']);
end

z=real(O);
end