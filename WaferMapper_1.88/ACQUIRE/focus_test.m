%% MAPFoSt test file
global GuiGlobalsStruct;
sm = GuiGlobalsStruct.MyCZEMAPIClass; %To shorten calls to global API variables in this function
clc;
tic;
single=1;
for i=1
% hardcoded aberrations
if single
    A=0;
    T1=15;
    T2=-15;
else
    A=[i 1 0]; %[um % %]
    T1=[15 0 0]; %defocus in [um]
    T2=[-15 0 0];
end



FileName = 'E:\PNI-Images\Eli\testimages\temFoc.tif';
ImageHeightInPixels = 1024;
ImageWidthInPixels = 1024;
DwellTimeInMicroseconds = 2;
PixSize=8; %nm/pixel
FOV=PixSize*ImageHeightInPixels/1000;%um
Acc=5;
NA= 0.5596*height / ((Acc*1000)^0.5);
%NA= sqrt(40)*0.752 / (PixSize/1000* (Acc*1000)^0.5);
%pause duration after changing WD
frametime=sm.Get_ReturnTypeSingle('AP_FRAME_TIME')/1000;
frametime=1;
%%
%Reset initial WD using AFStartingWDd from Montage Parameters -- should not be needed
% sm.Set_PassedTypeSingle('AP_WD',GuiGlobalsStruct.MontageParameters.AFStartingWD);

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
    %save(['F:\' mat2str(round(A)) mat2str(round(T1)) mat2str(round(T2)) 'PixSize' num2str(PixSize)],'A','T1','T2','I1','I2','Anom','T1nom','T2nom','FOV','StartingWorkingDistance','StartingStigX','StartingStigY');
end

% out=[I1 I2];
% imshow(out);
end

%% MAPFoSt

sigma =mean([std(double(I1(:))), std(double(I2(:)))]); %sigma for real space
[Kx, Ky]=meshgrid((mod(0.5+[0:ImageWidthInPixels-1]/ImageWidthInPixels,1)-0.5)*(6.28/FOV),(mod(0.5+[0:ImageHeightInPixels-1]/ImageHeightInPixels,1)-0.5)*(6.28/FOV));
fI1=fft2(double(I1)); %image should have dimension 2^n for faster FFT
fI2=fft2(double(I2));
if single
    init=2;
else
    init=[0 0 0];
end


p.length=20;
p.method='BFGS';
p.verbosity=1;
p.MFEPLS = 30;   % Max Func Evals Per Line Search
p.MSR = 100;                % Max Slope Ratio default
O=minimize(init,@MAP,p,fI1,fI2,T1,T2,NA,sigma,Kx,Ky,single);
disp(MAP(A,fI1,fI2,T1,T2,NA,sigma,Kx,Ky,single)-MAP(O,fI1,fI2,T1,T2,NA,sigma,Kx,Ky,single));

% if single
%     % plot 1-D MAP
%     fout=[];
%     dfout=[];
%     temp=-20:20;
%     for j=temp
%         [f,df]=MAP(j,fI1,fI2,T1,T2,NA,sigma,Kx,Ky,single);
%         fout=[fout f ];
%         dfout=[dfout df];
%     end
%     h=figure;
%     plot(temp,fout);
%     title(['experimental -ln(P(A)) vs A for A=' num2str(A)]);
% end
%%
if max(abs(O))<80
    if single
        sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance-10^-6*real(O)); %change WD to testing defocus
        disp(['final WD: ' num2str(10^6*sm.Get_ReturnTypeSingle('AP_WD')) 'um']);
    else
        sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance+10^-6*real(O(1)));
        sm.Set_PassedTypeSingle('AP_STIG_X',ResultStigX+10^-6*real(O(2)));
        sm.Set_PassedTypeSingle('AP_STIG_Y',ResultStigY+10^-6*real(O(3)));
        disp(['final WD: ' num2str(10^6*sm.Get_ReturnTypeSingle('AP_WD')) 'um']);
        disp(['final StigX: ' num2str(sm.Get_ReturnTypeSingle('AP_STIG_X')) '%']);
        disp(['final StigY: ' num2str(sm.Get_ReturnTypeSingle('AP_STIG_Y')) '%']);
    end

end


toc;


%clear I1 I2



%% function test
FileName = 'E:\PNI-Images\Eli\testimages\temFoc.tif';
ImageHeightInPixels = 128;
ImageWidthInPixels = 128;
DwellTimeInMicroseconds = 2;
PixSize=64; %nm/pixel
FOV=PixSize*ImageHeightInPixels/1000; %um

sm.Execute('CMD_AUTO_FOCUS_FINE');
pause(0.5);
while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
    pause(0.02);
end
focusWD1=sm.Get_ReturnTypeSingle('AP_WD');
sm.Execute('CMD_AUTO_FOCUS_FINE');
pause(0.5);
while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
    pause(0.02);
end
focusWD2=sm.Get_ReturnTypeSingle('AP_WD');
sm.Execute('CMD_AUTO_FOCUS_FINE');
pause(0.5);
while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
    pause(0.02);
end
focusWD3=sm.Get_ReturnTypeSingle('AP_WD');
focusWD=median([focusWD1 focusWD2 focusWD3]);

A=randi([-20 20],1,100); %generate sequence of initial aberrations
out=zeros(4,100);
for i=1:50

sm.Set_PassedTypeSingle('AP_WD',focusWD+10^-6*A(i));
pause(0.5);

tic;
[z,finalWD,I1,I2]=MAPFoSt(ImageHeightInPixels,ImageWidthInPixels,DwellTimeInMicroseconds,FileName,FOV,3);
time=toc;
bright=sm.Get_ReturnTypeString('AP_BRIGHTNESS');
contrast=sm.Get_ReturnTypeString('AP_CONTRAST');
stigx=sm.Get_ReturnTypeString('AP_STIG_X');
stigy=sm.Get_ReturnTypeString('AP_STIG_Y');
out(1,i)=10^6*focusWD;
out(2,i)=10^6*finalWD;
out(3,i)=A(i);
out(4,i)=time;
save(['F:\' 'Precision test (' mat2str(i) ')PixSize' num2str(PixSize)],'I1','I2','bright','contrast','stigx','stigy');
end

out2=zeros(4,100);
% for i=1:100
% 
% sm.Set_PassedTypeSingle('AP_WD',focusWD+10^-6*A(i));
% pause(0.5);
% 
% tic;
% sm.Execute('CMD_AUTO_FOCUS_FINE');
% pause(0.5);
% while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
%     pause(0.02);
% end
% finalWD=sm.Get_ReturnTypeSingle('AP_WD');
% time=toc;
% 
% out2(1,i)=10^6*focusWD;
% out2(2,i)=10^6*finalWD;
% out2(3,i)=A(i);
% out2(4,i)=time;
% end

save(['F:\' 'Precision test PixSize' num2str(PixSize)],'out','out2','FOV','A');
%% time built in zeiss AF
tic;
sm.Execute('CMD_AUTO_FOCUS_FINE');

pause(0.5);
disp('Auto Focusing...');
while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
    pause(0.02);
end
toc;