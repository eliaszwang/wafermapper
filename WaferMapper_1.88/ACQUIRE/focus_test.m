global GuiGlobalsStruct;
sm = GuiGlobalsStruct.MyCZEMAPIClass; %To shorten calls to global API variables in this function
clc;
tic;
single=1;
for i=1:10
% hardcoded aberrations
if single
    A=i;
    T1=7;
    T2=-7;
else
    A=[i i/2 i/3];
    T1=[15 0 0]; %defocus in [um]
    T2=[-15 0 0];
end



FileName = 'E:\PNI-Images\Eli\testimages\temFoc.tif';
FileName2 = 'E:\PNI-Images\Eli\testimages\temFoc2.tif';
ImageHeightInPixels = 1024;
ImageWidthInPixels = 1024;
DwellTimeInMicroseconds = 2;
PixSize=8; %nm/pixel
FOV=PixSize*ImageHeightInPixels/1000;%um
Acc=5;
NA= sqrt(40)*0.752 / (PixSize* (Acc*1000)^0.5);
%pause duration after changing WD
frametime=sm.Get_ReturnTypeSingle('AP_FRAME_TIME')/1000;
%frametime=0.1;
%%
%autofocus (AS-AF) at beginning
if single
    sm.Execute('CMD_AUTO_FOCUS_FINE');
    pause(0.5);
    while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
        pause(0.02);
    end
else
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
end
StartingWorkingDistance=sm.Get_ReturnTypeSingle('AP_WD');
StartingStigX = sm.Get_ReturnTypeSingle('AP_STIG_X');
StartingStigY = sm.Get_ReturnTypeSingle('AP_STIG_Y');
disp(['autofocus WD: ' num2str(10^6*StartingWorkingDistance) 'um']);
pause(frametime);
%%
%set initial aberration

if single
    sm.Set_PassedTypeSingle('AP_WD',StartingWorkingDistance+10^-6*A);
else
    sm.Set_PassedTypeSingle('AP_WD',StartingWorkingDistance+10^-6*A(1));
    sm.Set_PassedTypeSingle('AP_STIG_X',StartingStigX+10^-6*A(2));
    sm.Set_PassedTypeSingle('AP_STIG_Y',StartingStigY+10^-6*A(3));
end
pause(frametime);

%begin MAPFoSt
disp('Beginning MAPFoSt...');
CurrentWorkingDistance = sm.Get_ReturnTypeSingle('AP_WD');
CurrentStigX = sm.Get_ReturnTypeSingle('AP_STIG_X');
CurrentStigY = sm.Get_ReturnTypeSingle('AP_STIG_Y');
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
    sm.Set_PassedTypeSingle('AP_STIG_X',CurrentStigX+10^-6*T1(2));
    sm.Set_PassedTypeSingle('AP_STIG_Y',CurrentStigY+10^-6*T1(3));
end
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
    sm.Set_PassedTypeSingle('AP_STIG_X',CurrentStigX+10^-6*T2(2));
    sm.Set_PassedTypeSingle('AP_STIG_Y',CurrentStigY+10^-6*T2(3));
end
pause(frametime);
sm.Set_PassedTypeSingle('AP_SCANROTATION',0);
T2WD=sm.Get_ReturnTypeSingle('AP_WD');
T2StigX = sm.Get_ReturnTypeSingle('AP_STIG_X');
T2StigY = sm.Get_ReturnTypeSingle('AP_STIG_Y');
disp(['second image WD: ' num2str(10^6*T2WD) 'um']);
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

%get actual aberrations
if single
    Anom=A;
    A=10^6*(CurrentWorkingDistance-StartingWorkingDistance);
    T1nom=T1;
    T1=10^6*(T1WD-CurrentWorkingDistance);
    T2nom=T2;
    T2=10^6*(T2WD-CurrentWorkingDistance);
    %save(['F:\' mat2str(round([A T1 T2])) 'PixSize' num2str(PixSize)],'A','T1','T2','I1','I2','Anom','T1nom','T2nom','FOV');
else
    Anom=A;
    A(1)=10^6*(CurrentWorkingDistance-StartingWorkingDistance);
    A(2)=10^6*(CurrentStigX-StartingStigX);
    A(3)=10^6*(CurrentStigY-StartingStigY);
    T1nom=T1;
    T1(1)=10^6*(T1WD-CurrentWorkingDistance);
    T1(2)=10^6*(T1StigX-CurrentStigX);
    T1(3)=10^6*(T1StigY-CurrentStigY);
    T2nom=T2;
    T2(1)=10^6*(T2WD-CurrentWorkingDistance);
    T2(2)=10^6*(T2StigX-CurrentStigX);
    T2(3)=10^6*(T2StigY-CurrentStigY);
    %save(['F:\' mat2str(round(A)) mat2str(round(T1)) mat2str(round(T2)) 'PixSize' num2str(PixSize)],'A','T1','T2','I1','I2','Anom','T1nom','T2nom','FOV');
end

out=[I1 I2];
%imshow(out);
end

%% MAPFoSt
% 
% sigma =mean([std(double(I1(:))), std(double(I2(:)))]); %sigma for real space
% [Kx, Ky]=meshgrid((circshift([0:width-1]/width,width/2,2)-0.5)*(6.28/FOV),(circshift([0:width-1]/width,width/2,2)-0.5)*(6.28/FOV)); %units are rad/um?
% fI1=fft2(I1); %image should have dimension 2^n for faster FFT
% fI2=fft2(I2);
% if single
%     init=2;
% else
%     init=[0 0 0];
% end
% 
% 
% p.length=20;
% p.method='BFGS';
% p.verbosity=1;
% p.MFEPLS = 30;   % Max Func Evals Per Line Search
% p.MSR = 100;                % Max Slope Ratio default
% O=minimize(init,@MAP,p,fI1,fI2,T1,T2,NA,sigma,Kx,Ky,single);
% 
% %%
% if max(abs(O))<80
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
% 

toc;

%clear I1 I2

