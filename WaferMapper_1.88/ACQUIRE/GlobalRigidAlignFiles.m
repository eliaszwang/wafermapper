function transforms=GlobalRigidAlignFiles(Files,MatFiles)

%% calculate points of interest on each section
Options.init_sample=2;
Options.octaves=3;
Options.tresh=.003;
Options.centerfrac=.25;
[points,PixelRegion]=getSURFpointsfromFiles(Files,Options);
Z=length(points);


%% find relative similarity transforms between sections
% currently using ransac on top 50 matches.. see inside for details  
AlignmentOptions.verbose=0;
AlignmentOptions.PixelRegion=PixelRegion;
AlignmentOptions.ref_sect=25;
AlignmentOptions.dist_thresh=.02;
AlignmentOptions.min_inliers=20;
if exist('MatFiles','var')
    [transforms,num_inliers]=MetaAlign_SmallestBridge(points,AlignmentOptions,Files,MatFiles);
    %[transforms,numinliers,Pos1s,Pos2s,Inliers]=MetaAlignment(points,AlignmentOptions,Files,MatFiles);    
else
    [transforms,numinliers,Pos1s,Pos2s,Inliers]=MetaAlignment(points,AlignmentOptions,Files);   
end



% 
% figure(50);clf;
% imagesc(numinliers);colorbar;
% %% extract the relative parameters from the set of transforms
% [relscales,relangles,reldx,reldy,goodmatrix]=extract_relative_rigid(transforms,numinliers);
% 
% %% use gradient descent to optimize the angles from the relative angles
% [angles,angle_err,angle_Esave]=gradient_optimize_linear(relangles,2200,.01,goodmatrix);
% 
% %% use gradient descent to optimize the xpositions from the relative x shifts
% [xpos,xpos_err,xpos_Esave]=gradient_optimize_linear(reldx,2200,.01,goodmatrix);
% 
% %% use gradient descent to optimize the ypositions from the relative y shifts
% [ypos,ypos_err,ypos_Esave]=gradient_optimize_linear(reldy,2200,.01,goodmatrix);
% 
% avg_inliers=sum(numinliers,1)./sum(goodmatrix,1);
% avg_inliers(isnan(avg_inliers))=0;

   