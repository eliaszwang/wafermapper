function [z,a_on,a_diag]=MAPFoSt()
% Elias Wang's implementation of MAPFoSt
% algorithm adapted from Jonas Binding (Low-Dosage Maximum-A-Posteriori Focusing and Stigmation)

global GuiGlobalsStruct;

sm = GuiGlobalsStruct.MyCZEMAPIClass; %To shorten calls to global API variables in this function
targetFocus = GuiGlobalsStruct.MontageParameters.IsTargetFocus;
WDResetThreshold = GuiGlobalsStruct.MontageParameters.WDResetThreshold;


if ~exist('focOptions','var')
   focOptions.IsDoQualCheck = 0;
   focOptions.QualityThreshold = 0;
end

ImageHeightInPixels = 100; %currently MAPFoSt only compatible with square images
ImageWidthInPixels = 100;
DwellTimeInMicroseconds = 0.2;
FOV = GuiGlobalsStruct.MontageParameters.TileFOV_microns;


AFscanRate = GuiGlobalsStruct.MontageParameters.AutofunctionScanrate;
AFImageStore= GuiGlobalsStruct.MontageParameters.AutoFunctionImageStore;

if exist(GuiGlobalsStruct.TempImagesDirectory,'dir')
    FileName = [GuiGlobalsStruct.TempImagesDirectory '\tempFoc.tif'];
else
    FileName = 'C:\temp\temFoc.tif';
end


s = FOV/ImageHeightInPixels/1000000; %scale meters per pixel


StartingMagForAF =GuiGlobalsStruct.MontageParameters.AutoFocusStartMag;
StartingMagForAS = GuiGlobalsStruct.MontageParameters.AutoFocusStartMag;
startScanRot = sm.Get_ReturnTypeSingle('AP_SCANROTATION');

%Reset initial WD using AFStartingWDd from Montage Parameters
sm.Set_PassedTypeSingle('AP_WD',GuiGlobalsStruct.MontageParameters.AFStartingWD);
CurrentWorkingDistance = sm.Get_ReturnTypeSingle('AP_WD');

%Reset initial stig values using StartingStigX and StartingStigY from
%Montage Parameters
startStigX = GuiGlobalsStruct.MontageParameters.StartingStigX;
startStigY = GuiGlobalsStruct.MontageParameters.StartingStigY;

stage_x = sm.Get_ReturnTypeSingle('AP_STAGE_AT_X');
stage_y = sm.Get_ReturnTypeSingle('AP_STAGE_AT_Y');
stage_z = sm.Get_ReturnTypeSingle('AP_STAGE_AT_Z');
stage_t = sm.Get_ReturnTypeSingle('AP_STAGE_AT_T');
stage_r = sm.Get_ReturnTypeSingle('AP_STAGE_AT_R');
stage_m = sm.Get_ReturnTypeSingle('AP_STAGE_AT_M');

%%Take first image

sm.Set_PassedTypeSingle('AP_SCANROTATION',0);

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

%%Take second image

sm.Set_PassedTypeSingle('AP_SCANROTATION',0);

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

%%Calculate fft of two images
fI1=fft2(I1);
fI2=fft2(I2);

%%Find MAP
A=minimize([0 0 0],@MAP,5,fI1,fI2,T1,T2,FOV,ImageHeightInPixels,ImageWidthInPixels);


%set new WD/Stig from algorithm
sm.Set_PassedTypeSingle('AP_WD',A(1));
sm.Set_PassedTypeSingle('AP_STIG_X',A(2));
sm.Set_PassedTypeSingle('AP_STIG_Y',A(3));