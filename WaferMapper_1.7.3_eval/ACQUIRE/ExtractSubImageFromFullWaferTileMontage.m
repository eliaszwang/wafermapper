function [ SubImage ] = ExtractSubImageFromFullWaferTileMontage(Mouse_r_IndexInFullMap, Mouse_c_IndexInFullMap, WidthSubImage, HeightSubImage)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global GuiGlobalsStruct;

IsDisplayForDebugging = false;


Mouse_TileR = 1+floor(Mouse_r_IndexInFullMap/GuiGlobalsStruct.FullMapData.ImageHeightInPixels);
Mouse_TileC = 1+floor(Mouse_c_IndexInFullMap/GuiGlobalsStruct.FullMapData.ImageWidthInPixels);

%determine which tile R,C the user clicked in
AccumulatorImage3x3 = Read3x3TileImages(Mouse_TileR, Mouse_TileC);

r_in3x3 = Mouse_r_IndexInFullMap - (Mouse_TileR-2)*GuiGlobalsStruct.FullMapData.ImageHeightInPixels;
c_in3x3 = Mouse_c_IndexInFullMap - (Mouse_TileC-2)*GuiGlobalsStruct.FullMapData.ImageHeightInPixels;

min_r = floor(r_in3x3-HeightSubImage/2);
min_c = floor(c_in3x3-WidthSubImage/2);
SubImage = AccumulatorImage3x3(min_r:(min_r+HeightSubImage), min_c:(min_c+WidthSubImage));

if IsDisplayForDebugging
    h_3x3fig = figure();
    subplot(1,2,1);
    imshow(AccumulatorImage3x3,[0,255]);
    
    ColorArray = [1 0 0];
    h1 = line([c_in3x3-20, c_in3x3+20],[r_in3x3, r_in3x3]);
    h2 = line([c_in3x3, c_in3x3],[r_in3x3-20, r_in3x3+20]);
    set(h1,'Color',ColorArray);
    set(h2,'Color',ColorArray);
    
    h_top = line([c_in3x3-WidthSubImage/2, c_in3x3+WidthSubImage/2],[r_in3x3-HeightSubImage/2, r_in3x3-HeightSubImage/2]);
    set(h_top,'Color',ColorArray);
    h_bottom = line([c_in3x3-WidthSubImage/2, c_in3x3+WidthSubImage/2],[r_in3x3+HeightSubImage/2, r_in3x3+HeightSubImage/2]);
    set(h_bottom,'Color',ColorArray);
    h_left = line([c_in3x3-WidthSubImage/2, c_in3x3-WidthSubImage/2],[r_in3x3-HeightSubImage/2, r_in3x3+HeightSubImage/2]);
    set(h_left,'Color',ColorArray);
    h_right = line([c_in3x3+WidthSubImage/2, c_in3x3+WidthSubImage/2],[r_in3x3-HeightSubImage/2, r_in3x3+HeightSubImage/2]);
    set(h_right,'Color',ColorArray);
    
    set(h_3x3fig, 'Position', get(0,'Screensize')); %This maximizes the window for best resolution
    
    
    subplot(1,2,2);
    imshow(SubImage,[0,255]);
end

end

