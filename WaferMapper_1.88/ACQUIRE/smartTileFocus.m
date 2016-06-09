function[focusPosition] = smartTileFocus(StartingMagForAF, IsPerformAutoStig, StartingMagForAS, focOptions)

%% focus position is returned as XY
%% Retrieve necessary info
global GuiGlobalsStruct;


sm = GuiGlobalsStruct.MyCZEMAPIClass; %To shorten calls to global API variables in this function
targetFocus = GuiGlobalsStruct.MontageParameters.IsTargetFocus;
WDResetThreshold = GuiGlobalsStruct.MontageParameters.WDResetThreshold;


if ~exist('focOptions','var')
   focOptions.IsDoQualCheck = 0;
   focOptions.QualityThreshold = 0;
end

ImageHeightInPixels = 100;
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
GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('AP_WD',GuiGlobalsStruct.MontageParameters.AFStartingWD);

CurrentWorkingDistance = sm.Get_ReturnTypeSingle('AP_WD');

%Reset initial stig values using StartingStigX and StartingStigY from
%Montage Parameters
startStigX = GuiGlobalsStruct.MontageParameters.StartingStigX;
startStigY = GuiGlobalsStruct.MontageParameters.StartingStigY;

stage_x = sm.Get_ReturnTypeSingle('AP_STAGE_AT_X');
stage_y = sm.Get_ReturnTypeSingle('AP_STAGE_AT_Y');
stage_z = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STAGE_AT_Z');
stage_t = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STAGE_AT_T');
stage_r = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STAGE_AT_R');
stage_m = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STAGE_AT_M');

%%Take overview image

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
        I = imread(FileName);
    catch MyException
        IsReadOK = false;
        pause(0.1);
    end
end


%%find signal
if targetFocus

horiz = abs(I(:,1:end-1) - I(:,2:end));
vert = abs(I(1:end-1,:) - I(2:end,:));
difI = I;
difI(1:end-1,:) = difI(1:end-1,:)+vert;
difI(:,1:end-1) = difI(:,1:end-1)+horiz;

difI(difI>250) = 0;

focKern = ones(3,5);
focKern = focKern/sum(focKern(:));
fI= conv2(double(difI),focKern,'same');

[domeY domeX] = ind2sub(size(fI),1:(size(fI,1)*size(fI,2)));
domeDist = sqrt((domeY-mean(domeY)).^2 + (domeX-mean(domeX)).^2);
dome = fI *0; 
dome(:) = domeDist.^1.5;
dome = dome * max(fI(:))/max(dome(:)); 
fIdomed = fI-dome;

%image(fI),pause(.01)

[y x] = find(fIdomed == max(fIdomed(:)),1);

yshift = y - ImageHeightInPixels/2;
xshift = x - ImageHeightInPixels/2;

%% move stage to signal position
pause(.01)

GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('AP_Mag',150);
GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('DP_EXT_SCAN_CONTROL',0);
focusPosition = [stage_x - xshift * s stage_y + yshift * s];
sm.MoveStage(focusPosition(1),focusPosition(2),stage_z,stage_t,stage_r,stage_m);
while(strcmp(GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeString('DP_STAGE_IS'),'Busy'))
    pause(.01)
end

wmBackLash
pause(1)
end
%% Start Autofocus

IsNeedToReleaseFromFibics = 1;
if IsNeedToReleaseFromFibics
    %*** START: This sequence is designed to release the SEM from Fibics control
    sm.Execute('CMD_AUTO_FOCUS_FINE');
    pause(0.5);
    sm.Execute('CMD_ABORT_AUTO');
    while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
        pause(0.02);
    end
    sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance);
    
    pause(0.1);
    %*** END
end


sm.Set_PassedTypeSingle('AP_Mag',StartingMagForAF);
%%%%Turns off Ext Scan control here
%Temporary hard code settings
sm.Set_PassedTypeSingle('DP_AutoFunction_ScanRate',AFscanRate);
sm.Set_PassedTypeSingle('DP_IMAGE_STORE',AFImageStore);

%%%%
%sm.Set_PassedTypeSingle('DP_BEAM_BLANKED',0);
%%%%
%%LongPauseForStability
%sm.Execute('MCL_jmFocusWoble')%run for 10 sec
pause(3)
%sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance);


sm.Execute('CMD_AUTO_FOCUS_FINE');

pause(0.5);
disp('Auto Focusing...');
while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
    pause(0.02);
end
pause(0.01);
ResultWD =sm.Get_ReturnTypeSingle('AP_WD');
ResultWD1 = ResultWD*1000;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RS Add
ResultStigX =sm.Get_ReturnTypeSingle('AP_STIG_X');
ResultStigY =sm.Get_ReturnTypeSingle('AP_STIG_Y');

if IsPerformAutoStig
    %*** Auto stig
    %%%%%%%%%%%%%%%
    for repeatStig = 1:1
        sm.Set_PassedTypeSingle('AP_Mag',StartingMagForAS / repeatStig);
        
        AutoWorkingDistance =sm.Get_ReturnTypeSingle('AP_WD');
    if (abs(AutoWorkingDistance - CurrentWorkingDistance) > WDResetThreshold)
        %If the WD has gone far enough from CurrentWorkingDistance
        %(=AFStartingWD),reset it to CurrentWorkingDistance (i.e.
        %AFStartingWD, which is the WD used at the beginning of this
        %calculation)
        sm.Set_PassedTypeSingle('AP_WD', CurrentWorkingDistance);
        
        %Reset stig values with uesr-specified StartingStigX and
        %StartingStigY
        sm.Set_PassedTypeSingle('AP_STIG_X',startStigX);
        sm.Set_PassedTypeSingle('AP_STIG_Y',startStigY);
        pause(.1)
    end
        %Temporary hard code settings
        sm.Set_PassedTypeSingle('AP_Mag',5000 / repeatStig);
        sm.Set_PassedTypeSingle('DP_AutoFunction_ScanRate',AFscanRate);
        sm.Set_PassedTypeSingle('DP_IMAGE_STORE',AFImageStore);
        pause(0.01);
        
        
        %LogFile_WriteLine('Starting autostig');
        sm.Execute('CMD_AUTO_STIG');
        pause(0.5);
        
        %%%%%%%%%%%
        
        disp('Auto Stig...');
        while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
            pause(.1);
        end
        pause(0.1);
        ResultStigX1 = sm.Get_ReturnTypeSingle('AP_STIG_X');
        ResultStigY1 = sm.Get_ReturnTypeSingle('AP_STIG_Y');
        stig_difference_X = abs(ResultStigX1 - ResultStigX);
        stig_difference_Y = abs(ResultStigY1 - ResultStigY);
        if (stig_difference_X < GuiGlobalsStruct.MontageParameters.StigResetThreshold) && (stig_difference_Y < GuiGlobalsStruct.MontageParameters.StigResetThreshold)
            break % break out of repeating stig
        else
            'Repeating stig'
            %Reset Stig Values
            sm.Set_PassedTypeSingle('AP_STIG_X',startStigX);
            sm.Set_PassedTypeSingle('AP_STIG_Y',startStigY);
        end
        
    end %repeat Stig
    
    %*** Auto focusP
    %%%%%%%%%%%%%%%%%%%%%
    
    
    
    %%%%
    IsNeedToReleaseFromFibics = 1;
% if IsNeedToReleaseFromFibics
%     %*** START: This sequence is designed to release the SEM from Fibics control
%     sm.Execute('CMD_AUTO_FOCUS_FINE');
%     pause(0.5);
%     sm.Execute('CMD_ABORT_AUTO');
%     while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
%         pause(0.02);
%     end
%     sm.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance);
% 
%     pause(0.1);
%     %*** END
% end
%%%%

    sm.Set_PassedTypeSingle('AP_Mag',StartingMagForAF);
    pause(0.1);    
    sm.Set_PassedTypeSingle('DP_AutoFunction_ScanRate',AFscanRate);
    sm.Set_PassedTypeSingle('DP_IMAGE_STORE',AFImageStore);
    
    %%%%
    pause(.2)
    %%%%
    
    sm.Execute('CMD_AUTO_FOCUS_FINE');
    pause(0.5);
    disp('Auto Focusing...');
    while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
        pause(.1);
    end
    pause(0.1);
    
    
    ResultWD =sm.Get_ReturnTypeSingle('AP_WD');
    
    
    sm.Set_PassedTypeSingle('DP_AutoFunction_ScanRate',AFscanRate);
    sm.Set_PassedTypeSingle('DP_IMAGE_STORE',AFImageStore)
    
    
    AutoWorkingDistance = sm.Get_ReturnTypeSingle('AP_WD');
    if (abs(AutoWorkingDistance - CurrentWorkingDistance) > WDResetThreshold)
       
        
        %Reset to initial WD for calculation, CurrentWorkingDistance
        %(=AFStartingWD)
       %%%%Need to fix in Montage GUI
        sm.Set_PassedTypeSingle('AP_WD', GuiGlobalsStruct.CurrentWorkingDistance );
       %%%% 
        
        %%%% temporary until Montage GUI is working
         CurrentWorkingDistance =AFStartingWD;
          sm.Set_PassedTypeSingle('AP_WD', CurrentWorkingDistance );
         %%%% 
        %Reset to initial stig values
        sm.Set_PassedTypeSingle('AP_STIG_X',startStigX);
        sm.Set_PassedTypeSingle('AP_STIG_Y',startStigY);
        pause(.1)
    end
    
end %end Autofocus AutoStig autofocus

%% Return to original settings
if targetFocus
sm.Set_PassedTypeSingle('AP_SCANROTATION',startScanRot);
sm.MoveStage(stage_x ,stage_y ,stage_z,stage_t,stage_r,stage_m);
while(strcmp(GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeString('DP_STAGE_IS'),'Busy'))
    pause(.01)
end
wmBackLash
end

%% Refocus
if focOptions.IsDoQualCheck
    focOptions.IsDoQualCheck;
    offsetX = [0 50]/1000000;
    offsetY = [0 0 ]/1000000;
    StageX = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STAGE_AT_X');
    StageY = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STAGE_AT_Y');
    stage_z = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STAGE_AT_Z');
    stage_t = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STAGE_AT_T');
    stage_r = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STAGE_AT_R');
    stage_m = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STAGE_AT_M');
    
    focOptions.IsDoQualCheck = false;
    
    for o = 1:2
        [q] = takeFocusImage(focOptions);
         fprintf('Quality check registered %0.5g.\n',q.quality);
        %LogFile_WriteLine(sprintf('Quality check registered %0.5g.',q.quality));
        
     
        
        if q.quality <= focOptions.QualityThreshold %if bad
            disp('Image failed quality check')
           % LogFile_WriteLine(sprintf('!!!!! Autofocus failed to reach threshold %0.5g.',focOptions.QualityThreshold));
            
            %%Move over
            
            GuiGlobalsStruct.MyCZEMAPIClass.MoveStage(StageX + offsetX(o),StageY + offsetY(o),stage_z,stage_t,stage_r,stage_m);
            while(strcmp(GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeString('DP_STAGE_IS'),'Busy'))
                pause(.02)
            end
            wmBackLash
            
            % Reset WD to initial value again
            GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('AP_WD', CurrentWorkingDistance);

            GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('AP_STIG_X',startStigX);
            GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('AP_STIG_Y',startStigY);
            pause(.5)
            
            disp('try new focus')
            
            Perform_AF_or_AFASAF(StartingMagForAF, IsPerformAutoStig, StartingMagForAS, focOptions);
            
        else %if quality ok
            break
        end
    end %repeat refocus
    
    
    GuiGlobalsStruct.MyCZEMAPIClass.MoveStage(StageX,StageY,stage_z,stage_t,stage_r,stage_m);
    while(strcmp(GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeString('DP_STAGE_IS'),'Busy'))
        pause(.02)
    end
    wmBackLash
end





