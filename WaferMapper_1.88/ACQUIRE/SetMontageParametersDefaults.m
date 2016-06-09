function[] = SetMontageParametersDefaults()

global GuiGlobalsStruct
%% Set Default Parameters
defaultMontageParameters.TileFOV_microns = 40.96;
defaultMontageParameters.TileWidth_pixels = 4096;
defaultMontageParameters.TileDwellTime_microseconds = 1;

defaultMontageParameters.MontageNorthAngle = 0;
defaultMontageParameters.NumberOfTileRows = 3;
defaultMontageParameters.NumberOfTileCols = 3;
defaultMontageParameters.PercentTileOverlap = 6;
defaultMontageParameters.XOffsetFromAlignTargetMicrons = 0;
defaultMontageParameters.YOffsetFromAlignTargetMicrons = 0;
defaultMontageParameters.MicronsPerPixel = (defaultMontageParameters.TileFOV_microns/defaultMontageParameters.TileWidth_pixels);

defaultMontageParameters.AF_X_Offset_Microns = 0;
defaultMontageParameters.AF_Y_Offset_Microns = 0;

defaultMontageParameters.AutoFocusStartMag = 25000;
defaultMontageParameters.IsPerformQualityCheckOnEveryAF = false;
defaultMontageParameters.AFQualityThreshold = 3;
defaultMontageParameters.IsPerformQualCheckAfterEachImage = false;
defaultMontageParameters.ImageQualityThreshold = 3;

defaultMontageParameters.NoAuto = false;
defaultMontageParameters.IsTargetFocus =1;
defaultMontageParameters.IsSingle_AF_ForWholeMontage = false;
defaultMontageParameters.IsSingle_AFASAF_ForWholeMontage = true;
defaultMontageParameters.IsAFOnEveryTile = false;
defaultMontageParameters.IsAFASAFOnEveryTile = false;
defaultMontageParameters.IsPlaneFit = false;
defaultMontageParameters.IsXFit = false;
defaultMontageParameters.Is4square = false;
defaultMontageParameters.RowDistBetweenAFPointsMicrons = 50;
defaultMontageParameters.ColDistBetweenAFPointsMicrons = 50;
defaultMontageParameters.AutofunctionScanrate = 1;
defaultMontageParameters.AutoFunctionImageStore = 0;
defaultMontageParameters.IBSCContrast = 24;
defaultMontageParameters.IBSCBrightness = 45;
defaultMontageParameters.ImageContrast = 25;
defaultMontageParameters.ImageBrightness = 44;
defaultMontageParameters.retakeFocusType = 1;

defaultMontageParameters.IsAcquireOverviewImage = false;
defaultMontageParameters.MontageOverviewImageFOV_microns = 409.6;
defaultMontageParameters.MontageOverviewImageWidth_pixels = 4096;
defaultMontageParameters.MontageOverviewImageHeight_pixels = 4096;
defaultMontageParameters.MontageOverviewImageDwellTime_microseconds = 1;

defaultMontageParameters.AFStartingWD = .007; %In meters
% defaultMontageParameters.StartingStigX = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STIG_X'); %In percent
% defaultMontageParameters.StartingStigY = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STIG_Y'); %In percent
defaultMontageParameters.StartingStigX = 0.0; %In percent
defaultMontageParameters.StartingStigY = 0.0; %In percent

defaultMontageParameters.StigResetThreshold = 1.0; %In percent
defaultMontageParameters.WDResetThreshold = 0.0005; %m

defaultMontageParameters.AFTestImageFOV_microns = 40.96;
defaultMontageParameters.AFTestImageWidth_pixels = 4096;
defaultMontageParameters.AFTestImageDwellTime_microseconds = 1.0;
defaultMontageParameters.IsAFTestSameAsTileParameters = 0;

defaultFields = fields(defaultMontageParameters);

if isfield(GuiGlobalsStruct,'MontageParameters')  %update old paramets
    currentFields = fields(GuiGlobalsStruct.MontageParameters);

    for i = 1:length(defaultFields) 
        if ~sum(cell2mat(regexp(currentFields,defaultFields{i})))
           GuiGlobalsStruct.MontageParameters = ...
               setfield(GuiGlobalsStruct.MontageParameters, defaultFields{i},getfield(defaultMontageParameters,defaultFields{i})) ;
        end
    end
else
    GuiGlobalsStruct.MontageParameters = defaultMontageParameters;
end



