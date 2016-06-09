function [points1,points2]=RemovePointsFromOverlapZone(points1,points2,PixelRegion,MatFile1,MatFile2)
Rows=PixelRegion{1};
Cols=PixelRegion{2};
Info1=load(MatFile1,'Info');
Info2=load(MatFile2,'Info');
Info1=Info1.Info;
Info2=Info2.Info;

%buffersize=0; %microns %for buffer region that i don't think works or is
%necessary

%define a set of points which are the corners of the images in
%micron coordinates
%start with pixels, where 1,1 is at the upper left
Frame1x=[Cols(1) Cols(2) Cols(2) Cols(1)];
Frame1y=[Rows(1) Rows(1) Rows(2) Rows(2)];
%subtract off half the width in pixels so the center is in the
%middle
Frame1x=Frame1x-(Info1.ImageHeightInPixels/2);
Frame1y=Frame1y-(Info1.ImageHeightInPixels/2);
%Frame 2 is the same at this point in calculation
Frame2x=Frame1x;
Frame2y=Frame1y;

%now convert these pixels relative to the center coordinates
%into micrometers, where positive Y is down, and positive X is
%to the left
scale_micron_per_pixel=Info1.FOV_microns/Info1.ImageWidthInPixels;
Frame1x=-scale_micron_per_pixel*Frame1x+Info1.StageX_Meters*10^6;
Frame1y=scale_micron_per_pixel*Frame1y+Info1.StageY_Meters*10^6;
Frame2x=-scale_micron_per_pixel*Frame2x+Info2.StageX_Meters*10^6;
Frame2y=scale_micron_per_pixel*Frame2y+Info2.StageY_Meters*10^6;


%Since we can assume the rectangles are of equal height and
%width, then the corner which is within the other rectangle
%is the corner which makes up the corner of the rectnagle which
%defines the overlapping region.
F1CornerIsInside_F2=inpolygon(Frame1x,Frame1y,Frame2x,Frame2y);
F2CornerIsInside_F1=inpolygon(Frame2x,Frame2y,Frame1x,Frame1y);
if (sum(F1CornerIsInside_F2)>0)
    %then we have overlapping coordinates
    %pull out the coordinates of the corners which define the
    %overlapping rectangle
    F1CornerX=Frame1x(F1CornerIsInside_F2);
    F1CornerY=Frame1y(F1CornerIsInside_F2);
    F2CornerX=Frame2x(F2CornerIsInside_F1);
    F2CornerY=Frame2y(F2CornerIsInside_F1);
    %for buffer region that i don't think is working or necessary
%     if F1CornerX<F2CornerX
%         OverlapXv=[F1CornerX-buffersize F2CornerX+buffersize F2CornerX+buffersize F1CornerX-buffersize];
%     else
%         OverlapXv=[F1CornerX+buffersize F2CornerX-buffersize F2CornerX-buffersize F1CornerX+buffersize];
%     end
%     if F1CornerY<F2CornerY
%          OverlapYv=[F1CornerY-buffersize F1CornerY+buffersize F2CornerY+buffersize F2CornerY-buffersize];
%     else
%          OverlapYv=[F1CornerY+buffersize F1CornerY-buffersize F2CornerY-buffersize F2CornerY+buffersize];
%     end
    
    OverlapXv=[F1CornerX F2CornerX F2CornerX F1CornerX];
    OverlapYv=[F1CornerY F1CornerY F2CornerY F2CornerY];
    
    %pull out the coordinates of the points
    
    p1x=[points1(:).x];
    p1y=[points1(:).y];
    p2x=[points2(:).x];
    p2y=[points2(:).y];
    
    %convert them to positions in microns
    p1x=-scale_micron_per_pixel*p1x;
    p1y=scale_micron_per_pixel*p1y;
    p2x=-scale_micron_per_pixel*p2x;
    p2y=scale_micron_per_pixel*p2y;
    p1x=p1x+Info1.StageX_Meters*10^6;
    p1y=p1y+Info1.StageY_Meters*10^6;
    p2x=p2x+Info2.StageX_Meters*10^6;
    p2y=p2y+Info2.StageY_Meters*10^6;
    
    %now lets find which of these poitns are inside the overlapping rectangle
    
    [points_from_1_that_overlap]=inpolygon(p1x,p1y,OverlapXv,OverlapYv);
    [points_from_2_that_overlap]=inpolygon(p2x,p2y,OverlapXv,OverlapYv);
    
    if 1==1
        figure(53);
        clf;
        patch(Frame1x,Frame1y,'r');
        hold on;
        patch(Frame2x,Frame2y,'g');
        patch(OverlapXv,OverlapYv,'b');
        axis equal;
        scatter(p1x(~points_from_1_that_overlap),p1y(~points_from_1_that_overlap),'kx');
        scatter(p2x(~points_from_2_that_overlap),p2y(~points_from_2_that_overlap),'yx');
        
        set(gca,'Ydir','reverse');
        set(gca,'Xdir','reverse');
    end
    
    
    points1=points1(~points_from_1_that_overlap);
    points2=points2(~points_from_2_that_overlap);
    
end
