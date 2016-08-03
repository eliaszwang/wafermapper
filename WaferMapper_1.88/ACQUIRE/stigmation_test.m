%% MAPFoSt stigmation test
global GuiGlobalsStruct;
sm = GuiGlobalsStruct.MyCZEMAPIClass; %To shorten calls to global API variables in this function
clc;
tic;
single=0;


FileName = 'E:\PNI-Images\Eli\testimages\temFoc.tif';
ImageHeightInPixels = 1024;
ImageWidthInPixels = 1024;
DwellTimeInMicroseconds = 2;
PixSize=8; %nm/pixel
FOV=PixSize*ImageHeightInPixels/1000;%um
Acc=5;
%NA= 0.5596*ImageHeightInPixels / ((Acc*1000)^0.5);
%NA= sqrt(40)*0.752 / (PixSize/1000* (Acc*1000)^0.5);
%pause duration after changing WD
frametime=sm.Get_ReturnTypeSingle('AP_FRAME_TIME')/1000;
frametime=1;
%% defocusx
%autofocus (AS-AF) at beginning
% if single
%     sm.Execute('CMD_AUTO_FOCUS_FINE');
%     pause(0.5);
%     while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
%         pause(0.02);
%     end
% else
%     sm.Execute('CMD_AUTO_FOCUS_FINE');
%     pause(0.5);
%     while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
%         pause(0.02);
%     end
% %     sm.Execute('CMD_AUTO_STIG');
% %     pause(0.5);
% %     while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
% %         pause(.1);
% %     end
% %     pause(0.1);
% %     sm.Execute('CMD_AUTO_FOCUS_FINE');
% %     pause(0.5);
% %     while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
% %         pause(0.02);
% %     end
% end
StartingWorkingDistance=sm.Get_ReturnTypeSingle('AP_WD');
StartingStigX = sm.Get_ReturnTypeSingle('AP_STIG_X');
StartingStigY = sm.Get_ReturnTypeSingle('AP_STIG_Y');
disp(['autofocus WD: ' num2str(10^6*StartingWorkingDistance) 'um']);
pause(frametime);
for i=-20:20
% hardcoded aberrations
if single
    A=0;
    T1=15;
    T2=-15;
else
    A=[i 0 0]; %[um % %]
    T1=[15 0 0]; %defocus in [um]
    T2=[-15 0 0];
    T3=[7.5 0 0];
    T4=[-7.5 0 0];
end

%%
%set initial aberration

if single
    sm.Set_PassedTypeSingle('AP_WD',StartingWorkingDistance+10^-6*A);
else
    sm.Set_PassedTypeSingle('AP_WD',StartingWorkingDistance+10^-6*A(1));
    sm.Set_PassedTypeSingle('AP_STIG_X',StartingStigX+A(2));
    sm.Set_PassedTypeSingle('AP_STIG_Y',StartingStigY+A(3));
end
pause(frametime);

%begin MAPFoSt
disp('Beginning MAPFoSt...');
CurrentWorkingDistance = sm.Get_ReturnTypeSingle('AP_WD');
CurrentStigX = sm.Get_ReturnTypeSingle('AP_STIG_X');
CurrentStigY = sm.Get_ReturnTypeSingle('AP_STIG_Y');
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
if single
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T1);
else
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T1(1));
    sm.Set_PassedTypeSingle('AP_STIG_X',CurrentStigX+T1(2));
    sm.Set_PassedTypeSingle('AP_STIG_Y',CurrentStigY+T1(3));
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
    


%%Take second image
%implement errorcheck for test aberration
if single
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T2);
else
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T2(1));
    sm.Set_PassedTypeSingle('AP_STIG_X',CurrentStigX+T2(2));
    sm.Set_PassedTypeSingle('AP_STIG_Y',CurrentStigY+T2(3));
end
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


%%Take third image
%implement errorcheck for test aberration
if single
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T3);
else
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T3(1));
    sm.Set_PassedTypeSingle('AP_STIG_X',CurrentStigX+T3(2));
    sm.Set_PassedTypeSingle('AP_STIG_Y',CurrentStigY+T3(3));
end
pause(frametime);
sm.Set_PassedTypeSingle('AP_SCANROTATION',0);
T3WD=sm.Get_ReturnTypeSingle('AP_WD');
T3StigX = sm.Get_ReturnTypeSingle('AP_STIG_X');
T3StigY = sm.Get_ReturnTypeSingle('AP_STIG_Y');
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
        I3 = imread(FileName);
    catch MyException
        IsReadOK = false;
        pause(0.1);
    end
end

delete(FileName);


%%Take fourth image
%implement errorcheck for test aberration
if single
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T4);
else
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T4(1));
    sm.Set_PassedTypeSingle('AP_STIG_X',CurrentStigX+T4(2));
    sm.Set_PassedTypeSingle('AP_STIG_Y',CurrentStigY+T4(3));
end
pause(frametime);
sm.Set_PassedTypeSingle('AP_SCANROTATION',0);
T4WD=sm.Get_ReturnTypeSingle('AP_WD');
T4StigX = sm.Get_ReturnTypeSingle('AP_STIG_X');
T4StigY = sm.Get_ReturnTypeSingle('AP_STIG_Y');
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
        I4 = imread(FileName);
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
    A(2)=(CurrentStigX-StartingStigX);
    A(3)=(CurrentStigY-StartingStigY);
    T1nom=T1;
    T1(1)=10^6*(T1WD-CurrentWorkingDistance);
    T1(2)=(T1StigX-CurrentStigX);
    T1(3)=(T1StigY-CurrentStigY);
    T2nom=T2;
    T2(1)=10^6*(T2WD-CurrentWorkingDistance);
    T2(2)=(T2StigX-CurrentStigX);
    T2(3)=(T2StigY-CurrentStigY);
    save(['E:\PNI-Images\Eli\stigmation3\defocusx' mat2str(Anom) mat2str(T1nom) mat2str(T2nom) 'PixSize' num2str(PixSize)],'A','T1','T2','I1','I2','Anom','T1nom','T2nom','FOV','StartingWorkingDistance','StartingStigX','StartingStigY');
end

% out=[I1 I2];
% imshow(out);
end

%% defocusy
%autofocus (AS-AF) at beginning
% if single
%     sm.Execute('CMD_AUTO_FOCUS_FINE');
%     pause(0.5);
%     while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
%         pause(0.02);
%     end
% else
%     sm.Execute('CMD_AUTO_FOCUS_FINE');
%     pause(0.5);
%     while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
%         pause(0.02);
%     end
% %     sm.Execute('CMD_AUTO_STIG');
% %     pause(0.5);
% %     while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
% %         pause(.1);
% %     end
% %     pause(0.1);
% %     sm.Execute('CMD_AUTO_FOCUS_FINE');
% %     pause(0.5);
% %     while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
% %         pause(0.02);
% %     end
% end
% StartingWorkingDistance=sm.Get_ReturnTypeSingle('AP_WD');
% StartingStigX = sm.Get_ReturnTypeSingle('AP_STIG_X');
% StartingStigY = sm.Get_ReturnTypeSingle('AP_STIG_Y');
disp(['autofocus WD: ' num2str(10^6*StartingWorkingDistance) 'um']);
pause(frametime);
for i=-20:20
% hardcoded aberrations
if single
    A=0;
    T1=15;
    T2=-15;
else
    A=[i 1 1]; %[um % %]
    T1=[15 0 0]; %defocus in [um]
    T2=[-15 0 0];
end

%%
%set initial aberration

if single
    sm.Set_PassedTypeSingle('AP_WD',StartingWorkingDistance+10^-6*A);
else
    sm.Set_PassedTypeSingle('AP_WD',StartingWorkingDistance+10^-6*A(1));
    sm.Set_PassedTypeSingle('AP_STIG_X',StartingStigX+A(2));
    sm.Set_PassedTypeSingle('AP_STIG_Y',StartingStigY+A(3));
end
pause(frametime);

%begin MAPFoSt
disp('Beginning MAPFoSt...');
CurrentWorkingDistance = sm.Get_ReturnTypeSingle('AP_WD');
CurrentStigX = sm.Get_ReturnTypeSingle('AP_STIG_X');
CurrentStigY = sm.Get_ReturnTypeSingle('AP_STIG_Y');
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
if single
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T1);
else
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T1(1));
    sm.Set_PassedTypeSingle('AP_STIG_X',CurrentStigX+T1(2));
    sm.Set_PassedTypeSingle('AP_STIG_Y',CurrentStigY+T1(3));
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
    


%%Take second image
%implement errorcheck for test aberration
if single
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T2);
else
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T2(1));
    sm.Set_PassedTypeSingle('AP_STIG_X',CurrentStigX+T2(2));
    sm.Set_PassedTypeSingle('AP_STIG_Y',CurrentStigY+T2(3));
end
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

%%Take third image
%implement errorcheck for test aberration
if single
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T3);
else
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T3(1));
    sm.Set_PassedTypeSingle('AP_STIG_X',CurrentStigX+T3(2));
    sm.Set_PassedTypeSingle('AP_STIG_Y',CurrentStigY+T3(3));
end
pause(frametime);
sm.Set_PassedTypeSingle('AP_SCANROTATION',0);
T3WD=sm.Get_ReturnTypeSingle('AP_WD');
T3StigX = sm.Get_ReturnTypeSingle('AP_STIG_X');
T3StigY = sm.Get_ReturnTypeSingle('AP_STIG_Y');
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
        I3 = imread(FileName);
    catch MyException
        IsReadOK = false;
        pause(0.1);
    end
end

delete(FileName);


%%Take fourth image
%implement errorcheck for test aberration
if single
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T4);
else
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T4(1));
    sm.Set_PassedTypeSingle('AP_STIG_X',CurrentStigX+T4(2));
    sm.Set_PassedTypeSingle('AP_STIG_Y',CurrentStigY+T4(3));
end
pause(frametime);
sm.Set_PassedTypeSingle('AP_SCANROTATION',0);
T4WD=sm.Get_ReturnTypeSingle('AP_WD');
T4StigX = sm.Get_ReturnTypeSingle('AP_STIG_X');
T4StigY = sm.Get_ReturnTypeSingle('AP_STIG_Y');
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
        I4 = imread(FileName);
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
    A(2)=(CurrentStigX-StartingStigX);
    A(3)=(CurrentStigY-StartingStigY);
    T1nom=T1;
    T1(1)=10^6*(T1WD-CurrentWorkingDistance);
    T1(2)=(T1StigX-CurrentStigX);
    T1(3)=(T1StigY-CurrentStigY);
    T2nom=T2;
    T2(1)=10^6*(T2WD-CurrentWorkingDistance);
    T2(2)=(T2StigX-CurrentStigX);
    T2(3)=(T2StigY-CurrentStigY);
    save(['E:\PNI-Images\Eli\stigmation3\defocusy' mat2str(Anom) mat2str(T1nom) mat2str(T2nom) 'PixSize' num2str(PixSize)],'A','T1','T2','I1','I2','Anom','T1nom','T2nom','FOV','StartingWorkingDistance','StartingStigX','StartingStigY');
end

% out=[I1 I2];
% imshow(out);
end

%% stigx
%autofocus (AS-AF) at beginning
% if single
%     sm.Execute('CMD_AUTO_FOCUS_FINE');
%     pause(0.5);
%     while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
%         pause(0.02);
%     end
% else
%     sm.Execute('CMD_AUTO_FOCUS_FINE');
%     pause(0.5);
%     while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
%         pause(0.02);
%     end
% %     sm.Execute('CMD_AUTO_STIG');
% %     pause(0.5);
% %     while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
% %         pause(.1);
% %     end
% %     pause(0.1);
% %     sm.Execute('CMD_AUTO_FOCUS_FINE');
% %     pause(0.5);
% %     while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
% %         pause(0.02);
% %     end
% end
% StartingWorkingDistance=sm.Get_ReturnTypeSingle('AP_WD');
% StartingStigX = sm.Get_ReturnTypeSingle('AP_STIG_X');
% StartingStigY = sm.Get_ReturnTypeSingle('AP_STIG_Y');
disp(['autofocus WD: ' num2str(10^6*StartingWorkingDistance) 'um']);
pause(frametime);
for i=-2:0.1:2
% hardcoded aberrations
if single
    A=0;
    T1=15;
    T2=-15;
else
    A=[0 i 0]; %[um % %]
    T1=[15 0 0]; %defocus in [um]
    T2=[-15 0 0];
end

%%
%set initial aberration

if single
    sm.Set_PassedTypeSingle('AP_WD',StartingWorkingDistance+10^-6*A);
else
    sm.Set_PassedTypeSingle('AP_WD',StartingWorkingDistance+10^-6*A(1));
    sm.Set_PassedTypeSingle('AP_STIG_X',StartingStigX+A(2));
    sm.Set_PassedTypeSingle('AP_STIG_Y',StartingStigY+A(3));
end
pause(frametime);

%begin MAPFoSt
disp('Beginning MAPFoSt...');
CurrentWorkingDistance = sm.Get_ReturnTypeSingle('AP_WD');
CurrentStigX = sm.Get_ReturnTypeSingle('AP_STIG_X');
CurrentStigY = sm.Get_ReturnTypeSingle('AP_STIG_Y');
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
if single
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T1);
else
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T1(1));
    sm.Set_PassedTypeSingle('AP_STIG_X',CurrentStigX+T1(2));
    sm.Set_PassedTypeSingle('AP_STIG_Y',CurrentStigY+T1(3));
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
    


%%Take second image
%implement errorcheck for test aberration
if single
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T2);
else
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T2(1));
    sm.Set_PassedTypeSingle('AP_STIG_X',CurrentStigX+T2(2));
    sm.Set_PassedTypeSingle('AP_STIG_Y',CurrentStigY+T2(3));
end
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

%%Take third image
%implement errorcheck for test aberration
if single
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T3);
else
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T3(1));
    sm.Set_PassedTypeSingle('AP_STIG_X',CurrentStigX+T3(2));
    sm.Set_PassedTypeSingle('AP_STIG_Y',CurrentStigY+T3(3));
end
pause(frametime);
sm.Set_PassedTypeSingle('AP_SCANROTATION',0);
T3WD=sm.Get_ReturnTypeSingle('AP_WD');
T3StigX = sm.Get_ReturnTypeSingle('AP_STIG_X');
T3StigY = sm.Get_ReturnTypeSingle('AP_STIG_Y');
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
        I3 = imread(FileName);
    catch MyException
        IsReadOK = false;
        pause(0.1);
    end
end

delete(FileName);


%%Take fourth image
%implement errorcheck for test aberration
if single
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T4);
else
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T4(1));
    sm.Set_PassedTypeSingle('AP_STIG_X',CurrentStigX+T4(2));
    sm.Set_PassedTypeSingle('AP_STIG_Y',CurrentStigY+T4(3));
end
pause(frametime);
sm.Set_PassedTypeSingle('AP_SCANROTATION',0);
T4WD=sm.Get_ReturnTypeSingle('AP_WD');
T4StigX = sm.Get_ReturnTypeSingle('AP_STIG_X');
T4StigY = sm.Get_ReturnTypeSingle('AP_STIG_Y');
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
        I4 = imread(FileName);
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
    A(2)=(CurrentStigX-StartingStigX);
    A(3)=(CurrentStigY-StartingStigY);
    T1nom=T1;
    T1(1)=10^6*(T1WD-CurrentWorkingDistance);
    T1(2)=(T1StigX-CurrentStigX);
    T1(3)=(T1StigY-CurrentStigY);
    T2nom=T2;
    T2(1)=10^6*(T2WD-CurrentWorkingDistance);
    T2(2)=(T2StigX-CurrentStigX);
    T2(3)=(T2StigY-CurrentStigY);
    save(['E:\PNI-Images\Eli\stigmation3\stigx' mat2str((10*Anom)) mat2str(T1nom) mat2str(T2nom) 'PixSize' num2str(PixSize)],'A','T1','T2','I1','I2','Anom','T1nom','T2nom','FOV','StartingWorkingDistance','StartingStigX','StartingStigY');
end

% out=[I1 I2];
% imshow(out);
end


%% stigy
%autofocus (AS-AF) at beginning
% if single
%     sm.Execute('CMD_AUTO_FOCUS_FINE');
%     pause(0.5);
%     while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
%         pause(0.02);
%     end
% else
%     sm.Execute('CMD_AUTO_FOCUS_FINE');
%     pause(0.5);
%     while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
%         pause(0.02);
%     end
% %     sm.Execute('CMD_AUTO_STIG');
% %     pause(0.5);
% %     while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
% %         pause(.1);
% %     end
% %     pause(0.1);
% %     sm.Execute('CMD_AUTO_FOCUS_FINE');
% %     pause(0.5);
% %     while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
% %         pause(0.02);
% %     end
% end
% StartingWorkingDistance=sm.Get_ReturnTypeSingle('AP_WD');
% StartingStigX = sm.Get_ReturnTypeSingle('AP_STIG_X');
% StartingStigY = sm.Get_ReturnTypeSingle('AP_STIG_Y');
disp(['autofocus WD: ' num2str(10^6*StartingWorkingDistance) 'um']);
pause(frametime);

for i=-2:0.1:2
% hardcoded aberrations
if single
    A=0;
    T1=15;
    T2=-15;
else
    A=[0 0 i]; %[um % %]
    T1=[15 0 0]; %defocus in [um]
    T2=[-15 0 0];
end




%%
%set initial aberration

if single
    sm.Set_PassedTypeSingle('AP_WD',StartingWorkingDistance+10^-6*A);
else
    sm.Set_PassedTypeSingle('AP_WD',StartingWorkingDistance+10^-6*A(1));
    sm.Set_PassedTypeSingle('AP_STIG_X',StartingStigX+A(2));
    sm.Set_PassedTypeSingle('AP_STIG_Y',StartingStigY+A(3));
end
pause(frametime);

%begin MAPFoSt
disp('Beginning MAPFoSt...');
CurrentWorkingDistance = sm.Get_ReturnTypeSingle('AP_WD');
CurrentStigX = sm.Get_ReturnTypeSingle('AP_STIG_X');
CurrentStigY = sm.Get_ReturnTypeSingle('AP_STIG_Y');
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
if single
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T1);
else
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T1(1));
    sm.Set_PassedTypeSingle('AP_STIG_X',CurrentStigX+T1(2));
    sm.Set_PassedTypeSingle('AP_STIG_Y',CurrentStigY+T1(3));
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
    


%%Take second image
%implement errorcheck for test aberration
if single
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T2);
else
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T2(1));
    sm.Set_PassedTypeSingle('AP_STIG_X',CurrentStigX+T2(2));
    sm.Set_PassedTypeSingle('AP_STIG_Y',CurrentStigY+T2(3));
end
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

%%Take third image
%implement errorcheck for test aberration
if single
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T3);
else
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T3(1));
    sm.Set_PassedTypeSingle('AP_STIG_X',CurrentStigX+T3(2));
    sm.Set_PassedTypeSingle('AP_STIG_Y',CurrentStigY+T3(3));
end
pause(frametime);
sm.Set_PassedTypeSingle('AP_SCANROTATION',0);
T3WD=sm.Get_ReturnTypeSingle('AP_WD');
T3StigX = sm.Get_ReturnTypeSingle('AP_STIG_X');
T3StigY = sm.Get_ReturnTypeSingle('AP_STIG_Y');
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
        I3 = imread(FileName);
    catch MyException
        IsReadOK = false;
        pause(0.1);
    end
end

delete(FileName);


%%Take fourth image
%implement errorcheck for test aberration
if single
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T4);
else
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*T4(1));
    sm.Set_PassedTypeSingle('AP_STIG_X',CurrentStigX+T4(2));
    sm.Set_PassedTypeSingle('AP_STIG_Y',CurrentStigY+T4(3));
end
pause(frametime);
sm.Set_PassedTypeSingle('AP_SCANROTATION',0);
T4WD=sm.Get_ReturnTypeSingle('AP_WD');
T4StigX = sm.Get_ReturnTypeSingle('AP_STIG_X');
T4StigY = sm.Get_ReturnTypeSingle('AP_STIG_Y');
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
        I4 = imread(FileName);
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
    A(2)=(CurrentStigX-StartingStigX);
    A(3)=(CurrentStigY-StartingStigY);
    T1nom=T1;
    T1(1)=10^6*(T1WD-CurrentWorkingDistance);
    T1(2)=(T1StigX-CurrentStigX);
    T1(3)=(T1StigY-CurrentStigY);
    T2nom=T2;
    T2(1)=10^6*(T2WD-CurrentWorkingDistance);
    T2(2)=(T2StigX-CurrentStigX);
    T2(3)=(T2StigY-CurrentStigY);
    save(['E:\PNI-Images\Eli\stigmation3\stigy' mat2str((10*Anom)) mat2str(T1nom) mat2str(T2nom) 'PixSize' num2str(PixSize)],'A','T1','T2','I1','I2','Anom','T1nom','T2nom','FOV','StartingWorkingDistance','StartingStigX','StartingStigY');
end

% out=[I1 I2];
% imshow(out);
end