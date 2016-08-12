global GuiGlobalsStruct;
sm = GuiGlobalsStruct.MyCZEMAPIClass; %To shorten calls to global API variables in this function
clc;

filepath='E:\PNI-Images\Eli\testimages\';

%autofocus (AS-AF) at beginning

%     sm.Execute('CMD_AUTO_FOCUS_FINE');
%     pause(0.5);
%     while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
%         pause(0.02);
%     end

StartingWorkingDistance=sm.Get_ReturnTypeSingle('AP_WD');
StartingStigX = sm.Get_ReturnTypeSingle('AP_STIG_X');
StartingStigY = sm.Get_ReturnTypeSingle('AP_STIG_Y');
disp(['autofocus WD: ' num2str(10^6*StartingWorkingDistance) 'um']);
%pause duration after changing WD
%frametime=sm.Get_ReturnTypeSingle('AP_FRAME_TIME')/1000;
frametime=1;
pause(frametime);


for i=-25:25
tic;
% hardcoded aberrations

A=i;
T1=15;
T2=-15;



FileName = [filepath 'temFoc.tif'];
ImageHeightInPixels = 1024;
ImageWidthInPixels = 1024;
DwellTimeInMicroseconds = 2;
PixSize=8; %nm/pixel
FOV=PixSize*ImageHeightInPixels/1000;%um

%%
%set initial aberration
sm.Set_PassedTypeSingle('AP_WD',StartingWorkingDistance+10^-6*A);
pause(frametime);

%begin MAPFoSt
disp('Beginning MAPFoSt...');
CurrentWorkingDistance = sm.Get_ReturnTypeSingle('AP_WD');
CurrentStigX = sm.Get_ReturnTypeSingle('AP_STIG_X');
CurrentStigY = sm.Get_ReturnTypeSingle('AP_STIG_Y');
disp(['starting WD: ' num2str(10^6*CurrentWorkingDistance) 'um']);


%%Take first image
sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T1);
pause(frametime);
sm.Set_PassedTypeSingle('AP_SCANROTATION',0);
T1WD=sm.Get_ReturnTypeSingle('AP_WD');
T1StigX = sm.Get_ReturnTypeSingle('AP_STIG_X');
T1StigY = sm.Get_ReturnTypeSingle('AP_STIG_Y');
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


%%Take second image
sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T2);
pause(frametime);
sm.Set_PassedTypeSingle('AP_SCANROTATION',0);
T2WD=sm.Get_ReturnTypeSingle('AP_WD');
T2StigX = sm.Get_ReturnTypeSingle('AP_STIG_X');
T2StigY = sm.Get_ReturnTypeSingle('AP_STIG_Y');
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
toc;

%get actual aberrations

Anom=A;
A=10^6*(CurrentWorkingDistance-StartingWorkingDistance);
T1nom=T1;
T1=10^6*(T1WD-CurrentWorkingDistance);
T2nom=T2;
T2=10^6*(T2WD-CurrentWorkingDistance);
save([filepath 'MAPFoSt-test-images\test images 8_9_16\' mat2str([Anom T1nom T2nom]) 'PixSize' num2str(PixSize)],'A','T1','T2','I1','I2','Anom','T1nom','T2nom','FOV');


end