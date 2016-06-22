global GuiGlobalsStruct;
sm = GuiGlobalsStruct.MyCZEMAPIClass; %To shorten calls to global API variables in this function
clc;

single=1;
for i=1:10
% hardcoded aberrations
if single
    A=i;
    T1=7;
    T2=-7;
else
    A=[0 0 0];
    T1=[5 0 0]; %defocus in [um]
    T2=[-5 0 0];
end



FileName = 'E:\PNI-Images\Eli\testimages\temFoc.tif';
FileName2 = 'E:\PNI-Images\Eli\testimages\temFoc2.tif';
ImageHeightInPixels = 1024;
ImageWidthInPixels = 1024;
DwellTimeInMicroseconds = 2;
FOV=8.511;
Acc=5;
%pause duration after changing WD
frametime=sm.Get_ReturnTypeSingle('AP_FRAME_TIME')/1000;
%frametime=0.1;
%%
%autofocus (AS-AF) at beginning
sm.Execute('CMD_AUTO_FOCUS_FINE');
pause(0.5);
while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
    pause(0.02);
end
sm.Execute('CMD_AUTO_STIG');
pause(0.5);
while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
    pause(.1);
end
pause(0.1);
sm.Execute('CMD_AUTO_FOCUS_FINE');
pause(0.5);
while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
    pause(0.02);
end
CurrentWorkingDistance=sm.Get_ReturnTypeSingle('AP_WD');
ResultStigX = sm.Get_ReturnTypeSingle('AP_STIG_X');
ResultStigY = sm.Get_ReturnTypeSingle('AP_STIG_Y');
disp(['autofocus WD: ' num2str(10^6*CurrentWorkingDistance) 'um']);
pause(frametime);
%%
%set initial aberration

if single
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*A);
else
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*A(1));
    sm.Set_PassedTypeSingle('AP_STIG_X',ResultStigX+10^-6*A(2));
    sm.Set_PassedTypeSingle('AP_STIG_Y',ResultStigY+10^-6*A(3));
end
pause(frametime);

%begin MAPFoSt
disp('Beginning MAPFoSt...');
CurrentWorkingDistance = sm.Get_ReturnTypeSingle('AP_WD');
disp(['starting WD: ' num2str(10^6*CurrentWorkingDistance) 'um']);


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
   

%%Take first image
%implement errorcheck for test aberration
if single
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T1);
else
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T1(1));
    sm.Set_PassedTypeSingle('AP_STIG_X',ResultStigX+10^-6*T1(2));
    sm.Set_PassedTypeSingle('AP_STIG_Y',ResultStigY+10^-6*T1(3));
end
pause(frametime);
sm.Set_PassedTypeSingle('AP_SCANROTATION',0);
disp(['first image WD: ' num2str(10^6*sm.Get_ReturnTypeSingle('AP_WD')) 'um']);
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

%*** START: This sequence is designed to release the SEM from Fibics control
    sm.Execute('CMD_AUTO_FOCUS_FINE');
    pause(0.5);
    sm.Execute('CMD_ABORT_AUTO');
    while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
        pause(0.02);
    end
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance);
    
    pause(0.1);
    


%%Take second image
%implement errorcheck for test aberration
if single
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T2);
else
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T2(1));
    sm.Set_PassedTypeSingle('AP_STIG_X',ResultStigX+10^-6*T2(2));
    sm.Set_PassedTypeSingle('AP_STIG_Y',ResultStigY+10^-6*T2(3));
end
pause(frametime);
sm.Set_PassedTypeSingle('AP_SCANROTATION',0);
disp(['second image WD: ' num2str(10^6*sm.Get_ReturnTypeSingle('AP_WD')) 'um']);
sm.Fibics_WriteFOV(FOV);
%Wait for image to be acquired
sm.Fibics_AcquireImage(ImageWidthInPixels,ImageHeightInPixels,DwellTimeInMicroseconds,FileName2);
while(sm.Fibics_IsBusy)
    pause(.01); %1
end

%Wait for file to be written
IsReadOK = false;
while ~IsReadOK
    IsReadOK = true;
    try
        I2 = imread(FileName2);
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
    

save(['F:\' mat2str([A T1 T2])],'A','T1','T2','I1','I2');
out=[I1 I2];
%imshow(out);
end

%% MAPFoSt
% 
% 
% 
% if single
%     init=0;
%     MTF=@(Kx,Ky,A) exp(-0.125*(NA^2)*(Kx.^2+Ky.^2)*A^2);
% else
%     init=[0 0 0];
%     MTF=@(Kx,Ky,A) exp(-0.125*(NA^2)*(2*A(2)*(Kx.^2 - Ky.^2)*A(1) - A(3)*Kx.*Ky*A(1) + (Kx.^2+Ky.^2)*(A(2)^2 + A(3)^2 + A(1)^2)));
% end
% 
% 
% p.length=10;
% p.method='BFGS';
% p.verbosity=1;
% p.MFEPLS = 10;   % Max Func Evals Per Line Search
% p.MSR = 100;                % Max Slope Ratio default
% O=minimize(init,@MAP,p,I1,I2,T1,T2,FOV,Acc,single);
% 
% %%
% if max(real(O))<80
%     if single
%         sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance-10^-6*real(O)); %change WD to testing defocus
%         disp(['final WD: ' num2str(10^6*sm.Get_ReturnTypeSingle('AP_WD')) 'um']);
%     else
%         sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*real(O(1)));
%         sm.Set_PassedTypeSingle('AP_STIG_X',ResultStigX+10^-6*real(O(2)));
%         sm.Set_PassedTypeSingle('AP_STIG_Y',ResultStigY+10^-6*real(O(3)));
%         disp(['final WD: ' num2str(10^6*sm.Get_ReturnTypeSingle('AP_WD')) 'um']);
%         disp(['final StigX: ' num2str(sm.Get_ReturnTypeSingle('AP_STIG_X')) '%']);
%         disp(['final StigY: ' num2str(sm.Get_ReturnTypeSingle('AP_STIG_Y')) '%']);
%     end
% 
% end


clear I1 I2

