function [z,a_on,a_diag]=MAPFoSt(FOV,ImageHeightInPixels,ImageWidthInPixels,DwellTimeInMicroseconds,FileName,AccV)
% Elias Wang's implementation of MAPFoSt
% algorithm adapted from Jonas Binding (Low-Dosage Maximum-A-Posteriori Focusing and Stigmation)
% should replace "sm.Execute('CMD_AUTO_FOCUS_FINE')" calls
global GuiGlobalsStruct;
sm = GuiGlobalsStruct.MyCZEMAPIClass; %To shorten calls to global API variables in this function

%set AccV as optional for now
if ~exist('AccV','var')
    AccV=100; %just random number, need to be changed
end
    
% hardcoded test aberrations (defocus only)
T1=[15 0 0]; %defocus in [um]
T2=[-15 0 0];

%Reset initial WD using AFStartingWDd from Montage Parameters -- should not be needed
% sm.Set_PassedTypeSingle('AP_WD',GuiGlobalsStruct.MontageParameters.AFStartingWD);

%Get current values for WD and Stig
CurrentWorkingDistance = sm.Get_ReturnTypeSingle('AP_WD');
ResultStigX =sm.Get_ReturnTypeSingle('AP_STIG_X');
ResultStigY =sm.Get_ReturnTypeSingle('AP_STIG_Y');

%%Take first image
%implement errorcheck for test aberration
sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T1(1)); %change WD to testing defocus
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
%implement errorcheck for test aberration
sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T2(1)); %change WD to testing defocus
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

%% Find MAP
% initialize minimization at [0 0 0]
%set to 5 iterations for now
A=minimize([0 0 0],@MAP,5,I1,I2,T1,T2,FOV,AccV,1);%last parameter sets single(1)/multiple(0) aberration modes

% set new WD/Stig from algorithm
z=A(1);
a_on=A(2);
a_diag=A(3);
sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*A(1));
sm.Set_PassedTypeSingle('AP_STIG_X',ResultStigX+10^-6*A(2));
sm.Set_PassedTypeSingle('AP_STIG_Y',ResultStigY+10^-6*A(3));
end