global GuiGlobalsStruct;


sm = GuiGlobalsStruct.MyCZEMAPIClass; %To shorten calls to global API variables in this function


    FileName = 'E:\PNI-Images\Eli\testimages\temFoc.tif';

ImageHeightInPixels = 1024;
ImageWidthInPixels = 1024;
DwellTimeInMicroseconds = 2;
FOV=8.511;
sm.Execute('CMD_AUTO_FOCUS_FINE');

pause(0.5);
disp('Auto Focusing...');
while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
    pause(0.02);
end
pause(5);


CurrentWorkingDistance = sm.Get_ReturnTypeSingle('AP_WD');

% hardcoded test aberrations (defocus only)
T1=[0 0 0]; %defocus in [um]
T2=[0 0 0];

%Reset initial WD using AFStartingWDd from Montage Parameters -- should not be needed
% sm.Set_PassedTypeSingle('AP_WD',GuiGlobalsStruct.MontageParameters.AFStartingWD);
% CurrentWorkingDistance = sm.Get_ReturnTypeSingle('AP_WD');

%*** START: This sequence is designed to release the SEM from Fibics control
    sm.Execute('CMD_AUTO_FOCUS_FINE');
    pause(0.5);
    sm.Execute('CMD_ABORT_AUTO');
    while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
        pause(0.02);
    end
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance);
    
    pause(0.1);
    
    pause(3)
sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance);

%%Take first image
%implement errorcheck for test aberration
sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T1(1)); %change WD to testing defocus
sm.Set_PassedTypeSingle('AP_SCANROTATION',0);
sm.Get_ReturnTypeSingle('AP_WD')*1000
sm.Fibics_WriteFOV(FOV);
%Wait for image to be acquired
sm.Fibics_AcquireImage(ImageWidthInPixels,ImageHeightInPixels,DwellTimeInMicroseconds,FileName);
disp('getting first image');
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

%*** START: This sequence is designed to release the SEM from Fibics control
    sm.Execute('CMD_AUTO_FOCUS_FINE');
    pause(0.5);
    sm.Execute('CMD_ABORT_AUTO');
    while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
        pause(0.02);
    end
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance);
    
    pause(0.1);
    
    pause(3)
sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance);


%%Take second image
%implement errorcheck for test aberration
sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T2(1)); %change WD to testing defocus
sm.Set_PassedTypeSingle('AP_SCANROTATION',0);
sm.Get_ReturnTypeSingle('AP_WD')*1000
sm.Fibics_WriteFOV(FOV);
%Wait for image to be acquired
sm.Fibics_AcquireImage(ImageWidthInPixels,ImageHeightInPixels,DwellTimeInMicroseconds,FileName);
disp('getting second image');
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


%*** START: This sequence is designed to release the SEM from Fibics control
    sm.Execute('CMD_AUTO_FOCUS_FINE');
    pause(0.5);
    sm.Execute('CMD_ABORT_AUTO');
    while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
        pause(0.02);
    end
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance);
    
    pause(0.1);
    
    pause(3)
sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance);

out=[I1 I2];
imshow(out);

clear I1 I2;