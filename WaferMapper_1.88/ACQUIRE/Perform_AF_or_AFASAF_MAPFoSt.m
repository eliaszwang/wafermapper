function [] = Perform_AF_or_AFASAF_MAPFoSt(StartingMagForAF, IsPerformAutoStig, StartingMagForAS, focOptions, IsNeedToReleaseFromFibics)

global GuiGlobalsStruct;

AFscanRate = GuiGlobalsStruct.MontageParameters.AutofunctionScanrate;
AFImageStore = GuiGlobalsStruct.MontageParameters.AutoFunctionImageStore;

%GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('DP_IMAGE_STORE',0); %2

if ~exist('focOptions', 'var')
    focOptions.QualityThreshold = 0;
    focOptions.IsDoQualCheck = false;
end

% if focOptions.IsDoQualCheck % focus fast if you are going to check later
%     GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('DP_IMAGE_STORE',1);
% end

if ~exist('IsNeedToReleaseFromFibics', 'var')
    IsNeedToReleaseFromFibics = true;
end

%set up MAPFoSt variables
MAPFoSt_image_height=128; %test image height in pixels
MAPFoSt_image_width=128; %test image width in pixels
MAPFoSt_image_dwelltime_us=2; %test image dwelltime in microseconds
MAPFoSt_image_FOV=8.192; %test image FOV in microns
MAPFoSt_maxiter=3; %max number of iterations for MAPFoSt
MAPFoSt_fallback=1; %boolean for defaulting to Zeiss AF
MAPFoSt_verbosity=0; 

if exist(GuiGlobalsStruct.TempImagesDirectory,'dir')
    MAPFoSt_image_filename = [GuiGlobalsStruct.TempImagesDirectory '\tempFoc.tif']; %test image filename for saving
else
    MAPFoSt_image_filename = 'C:\temp\temFoc.tif';
end

CurrentWorkingDistance = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_WD');
% Change to the user-specified stigmation values StartingStigX,
% StartingStigY, instead of the current stigmation values
startStigX = GuiGlobalsStruct.MontageParameters.StartingStigX;
startStigY = GuiGlobalsStruct.MontageParameters.StartingStigY;



% if IsNeedToReleaseFromFibics
%     %*** START: This sequence is designed to release the SEM from Fibics control
%     GuiGlobalsStruct.MyCZEMAPIClass.Execute('CMD_AUTO_FOCUS_FINE');
%     pause(0.5);
%     GuiGlobalsStruct.MyCZEMAPIClass.Execute('CMD_ABORT_AUTO');
%     pause(0.5);
%    
%     GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('AP_WD',CurrentWorkingDistance);
%     pause(0.1);
%     %*** END
% end

%*** Auto focus
%%%%%%%%%%%%%%%%%%%%%%%%%%
fstart= tic;
%LogFile_WriteLine('   Starting first autofocus');



%Temporary hard code settings
GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('DP_AutoFunction_ScanRate',AFscanRate); %these settings might not be necessary, EW
GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('DP_IMAGE_STORE',AFImageStore);
GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('AP_Mag',StartingMagForAF);

%%%%
pause(1)
%%%%

disp('Auto Focusing...');
MAPFoSt(MAPFoSt_image_height,MAPFoSt_image_width,MAPFoSt_image_dwelltime_us,MAPFoSt_image_filename,MAPFoSt_image_FOV,MAPFoSt_maxiter,MAPFoSt_fallback,MAPFoSt_verbosity);

time_elapsed = toc(fstart);
ResultWD = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_WD');
ResultWD1 = ResultWD*1000; %In mm
%LogFile_WriteLine(sprintf('   Beginning WD =                %0.6g mm',CurrentWorkingDistance*1000));
%LogFile_WriteLine(sprintf('   First autofocus produced WD = %0.6g mm',ResultWD*1000));
%focus_difference = abs(ResultWD*1000 - CurrentWorkingDistance*1000)*1000;
% Commented out the line above because it is unused in this code
%LogFile_WriteLine(sprintf('   WD difference =               %0.5g um',focus_difference));
%LogFile_WriteLine(sprintf('   First autofocus time = %0.5g sec',time_elapsed));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RS Add
ResultStigX = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STIG_X');
ResultStigY = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STIG_Y');
%LogFile_WriteLine(sprintf('      Starting Stig X =     %0.5g and Y = %0.5g',ResultStigX, ResultStigY));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if IsPerformAutoStig
    %*** Auto stig
    %%%%%%%%%%%%%%%
    astig = tic;
    for repeatStig = 1:3
        GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('AP_Mag',StartingMagForAS / repeatStig);
        
        %***Fix these hard code settings***
        %Temporary hard code settings
        GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('AP_Mag',30000 / repeatStig);
        GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('DP_AutoFunction_ScanRate',AFscanRate);
        GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('DP_IMAGE_STORE',AFImageStore)
        pause(0.5);
        
        
        %LogFile_WriteLine('Starting autostig');
        GuiGlobalsStruct.MyCZEMAPIClass.Execute('CMD_AUTO_STIG');
        pause(0.5);
        
        
        %%%%%%%%%%%
        
        disp('Auto Stig...');
        while ~strcmp('Idle',GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeString('DP_AUTO_FUNCTION'))
            pause(1);
        end
        pause(0.5);
        ResultStigX1 = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STIG_X')
        ResultStigY1 = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STIG_Y')
        stig_time = toc(astig);
        stig_difference_X = abs(ResultStigX1 - ResultStigX);
        stig_difference_Y = abs(ResultStigY1 - ResultStigY);
        if (stig_difference_X < GuiGlobalsStruct.MontageParameters.StigResetThreshold) && (stig_difference_Y < GuiGlobalsStruct.MontageParameters.StigResetThreshold)
            break % break out of repeating stig
        else
            'Repeating stig'
            %Reset Stig Values
            GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('AP_STIG_X',startStigX);
            GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('AP_STIG_Y',startStigY);
        end
        
    end %repeat Stig
    %LogFile_WriteLine(sprintf('      AutoStig produced X = %0.5g  and Y = %0.5g',ResultStigX1, ResultStigY1));
    %LogFile_WriteLine(sprintf('      Diff Stig X =         %0.5g   and Y = %0.5g',stig_difference_X, stig_difference_Y));
    %LogFile_WriteLine(sprintf('      Stigmation Time = %0.5g sec',stig_time));
    
    
    %*** Auto focus
    %%%%%%%%%%%%%%%%%%%%%
    f2start = tic;
    %LogFile_WriteLine('        Starting second autofocus');
    GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('AP_Mag',StartingMagForAF); %these settings might not be necessary, EW
    %pause(0.5);
    
    
    GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('DP_AutoFunction_ScanRate',AFscanRate);
    GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('DP_IMAGE_STORE',AFImageStore)
    %%%%
        pause(1);
    %%%%
    
    disp('Auto Focusing...');
    MAPFoSt(MAPFoSt_image_height,MAPFoSt_image_width,MAPFoSt_image_dwelltime_us,MAPFoSt_image_filename,MAPFoSt_image_FOV,MAPFoSt_maxiter,MAPFoSt_fallback,MAPFoSt_verbosity);

    
    time_to_focus = toc(f2start);
    
    ResultWD = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_WD');
    
    %LogFile_WriteLine(sprintf('        Second autofocus produced WD = %0.5g mm',ResultWD*1000));
    focus_difference2 = abs(ResultWD - ResultWD1);
    %LogFile_WriteLine(sprintf('        Final WD difference =          %0.5g um',focus_difference2));
    %LogFile_WriteLine(sprintf('        Second autofocus Time = %0.5g sec',time_to_focus));
    ResultBrightness = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_Brightness');
    ResultContrast = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_Contrast');
    %LogFile_WriteLine(sprintf('           Brightness = %0.5g, Contrast = %0.5g',ResultBrightness, ResultContrast));
    
end %if stig
    
    
    %%%%%%%%%%%
    % ResultWD = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_WD');
    % LogFile_WriteLine(sprintf('Autofocus produced WD = %0.5g',ResultWD*1000));
    
    if focOptions.IsDoQualCheck
        
        StageX = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STAGE_AT_X');
        StageY = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STAGE_AT_Y');
        stage_z = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STAGE_AT_Z');
        stage_t = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STAGE_AT_T');
        stage_r = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STAGE_AT_R');
        stage_m = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STAGE_AT_M');
        
        focOptions.IsDoQualCheck = false;
        
        [q] = takeFocusImage(focOptions);
        %LogFile_WriteLine(sprintf('Quality check registered %0.5g.',q.quality))
        
        if q.quality <= focOptions.QualityThreshold %if bad
            'Image failed quality check'
            %LogFile_WriteLine(sprintf('!!!!! Autofocus failed to reach threshold %0.5g.',focOptions.QualityThreshold))
            fprintf('!!!!! Autofocus failed to reach threshold %0.5g.\n',focOptions.QualityThreshold)
            %%Move over
            
            GuiGlobalsStruct.MyCZEMAPIClass.MoveStage(StageX + offsetX(o),StageY + offsetY(o),stage_z,stage_t,stage_r,stage_m);
            while(strcmp(GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeString('DP_STAGE_IS'),'Busy'))
                pause(.02)
            end
            wmBackLash
            
            GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('AP_WD', GuiGlobalsStruct.MontageParameters.AFStartingWD);
            GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('AP_STIG_X',startStigX);
            GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('AP_STIG_Y',startStigY);
            pause(.2)
            
            'try new focus'
            
            %%Run smart focus to get focused on region with tissue
            smartTileFocus(StartingMagForAF, IsPerformAutoStig, StartingMagForAS, focOptions);
            
            
            
            GuiGlobalsStruct.MyCZEMAPIClass.MoveStage(StageX,StageY,stage_z,stage_t,stage_r,stage_m);
            while(strcmp(GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeString('DP_STAGE_IS'),'Busy'))
                pause(.02)
            end
            wmBackLash
        end
    end
    
        GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('DP_AutoFunction_ScanRate',AFscanRate);
        GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('DP_IMAGE_STORE',AFImageStore)
        
        
        AutoWorkingDistance = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_WD');
        if (abs(AutoWorkingDistance - CurrentWorkingDistance)) > (GuiGlobalsStruct.MontageParameters.WDResetThreshold) %NOTE: this was hard-coded at a value 3 times the normal threshold size (0.5mm).  Not sure why.    
            GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('AP_WD', GuiGlobalsStruct.MontageParameters.AFStartingWD);
            GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('AP_STIG_X',startStigX);
            GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('AP_STIG_Y',startStigY);
            pause(.1)
        end
        
    end
    
    
    
    
    %correct for stage move if done