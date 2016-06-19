global GuiGlobalsStruct;


sm = GuiGlobalsStruct.MyCZEMAPIClass;

  FileName = 'E:\PNI-Images\Eli\testimages\temFoc.tif';

ImageHeightInPixels = 1024;
ImageWidthInPixels = 1024;
DwellTimeInMicroseconds = 2;
CurrentWorkingDistance = sm.Get_ReturnTypeSingle('AP_WD');





[z,a_on,a_diag]=MAPFoSt(CurrentWorkingDistance,8.511,ImageHeightInPixels,ImageWidthInPixels,DwellTimeInMicroseconds,FileName,5);