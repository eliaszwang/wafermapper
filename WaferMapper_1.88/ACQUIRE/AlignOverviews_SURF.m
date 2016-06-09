function[] = AlignOverviews

global GuiGlobalsStruct

%pull out the list of files and Matfiles
%[Files,MatFiles,labels]=GetSortedImagesAndMatfiles(GuiGlobalsStruct.SectionOverviewsDirectory);




if GuiGlobalsStruct.SectionOverviewProcessingParameters.SURFall
        waferList = GuiGlobalsStruct.ListOfWaferNames;
       

else
    waferList{1} = GuiGlobalsStruct.WaferName;
end

secOverviewFiles = GetSortedImagesAndMatfiles_allWafers(waferList);
Files = secOverviewFiles.Files;
MatFiles = secOverviewFiles.MatFiles;
labels = secOverviewFiles.labels;

%[xpos,ypos,angles,avg_inliers]=GlobalRigidAlignFiles(Files,MatFiles);

%extract SURF points from each section
SurfOptions=GuiGlobalsStruct.SurfOptions;
RanSacOptions=GuiGlobalsStruct.RanSacOptions;
MetaAlignmentOptions=GuiGlobalsStruct.MetaAlignmentOptions;
[points,PixelRegion]=getSURFpointsfromFiles(Files,SurfOptions);

%now do something different depending on the alignment method
switch MetaAlignmentOptions.Method
    case{'GradDescent'} %not finished
        Options.verbose=0;
        
        Options.Nbest=RanSacOptions.NBest;
        Options.dist_thresh=RanSacOptions.dist_thresh;
        Options.det_thresh=RanSacOptions.MaxDetChange;
        
        Options.n_dense=MetaAlignmentOptions.LocalDepth;
        Options.long_length= MetaAlignmentOptions.LongRangeDelta;
        Options.n_long=MetaAlignmentOptions.LongRangeDepth;
        
        [transforms,num_inliers,Pos1s,Pos2s,Inliers]=MetaAlignment(points,Options,Files,MatFiles);
        [dtheta,dx,dy]=ExtractDthetaDxDx(transforms);
        
%         Z=size(num_inliers,1);
%         for i=1:Z
%             for j=1:Z
%                 
%             end
%         end
%         MetaAlignmentOptions.MinInliers
        
    case{'MatrixInversion'} %not finished
        Options.verbose=0;
        
        Options.Nbest=RanSacOptions.NBest;
        Options.dist_thresh=RanSacOptions.dist_thresh;
        Options.det_thresh=RanSacOptions.MaxDetChange;
        
        Options.n_dense=MetaAlignmentOptions.LocalDepth;
        Options.long_length= MetaAlignmentOptions.LongRangeDelta;
        Options.n_long=MetaAlignmentOptions.LongRangeDepth;
        Options.PixelRegion=PixelRegion;
        [transforms,num_inliers,Pos1s,Pos2s,Inliers]=MetaAlignment(points,Options,Files,MatFiles);
      
        [dTdpX,dTdpY]=MakeConstraintMatrix(num_inliers,Pos1s,Pos2s,15);
                
        lambda=100000;
        Z=size(num_inliers,1);
        C=zeros(3*Z,1);
        C(1:3:3*size(num_inliers,1))=lambda;
        M=(dTdpX + lambda*eye(3*Z,3*Z));
        
        P=C/lambda;
        
        
        delt=.00000001;
        for i=1:100
            dE=M*P-C;
            P=P-delt*dE;
            figure(3);
            clf;
            plot(P(1:3:end));
            pause(.1);
        end
        X=inv(M)*C;
        
    case{'ShortBridge'}
        Options.verbose=0;
        Options.Nbest=RanSacOptions.NBest;
        Options.dist_thresh=RanSacOptions.dist_thresh;
        Options.ref_sect=round(length(Files)/2);
        Options.det_thresh=RanSacOptions.MaxDetChange;
        Options.min_inliers=MetaAlignmentOptions.MinInliers;
        Options.max_dist=MetaAlignmentOptions.MaxDelta;
        Options.PixsizelRegion=PixelRegion;
        [transforms,num_inliers,notSectionFlag]=MetaAlign_SmallestBridge(points,Options,Files,MatFiles);
    otherwise
        disp('Error!');
end

thesize=matlabpool('size');
if thesize==0
    matlabpool OPEN;
end
Ms=cell(1,length(Files));

parfor i = 1:length(Files)
    disp(i);
 

    OverviewImageAlignedFileNameStr =[secOverviewFiles.alignedDir{i} 'SectionOverviewAligned_' num2str(labels(i)) '.tif'];

    OverviewImage=imread(Files{i});    
    M=ConvertCompositeAffineTransform(transforms{i});
    Ms{i}=M;
    tform=maketform('affine',M);
    OverviewImage_rotated_shifted=imtransform(OverviewImage,tform,'bicubic',...
        'UData',[-size(OverviewImage,2)/2 size(OverviewImage,2)/2],...
        'VData',[-size(OverviewImage,1)/2 size(OverviewImage,1)/2],...
        'XData',[-size(OverviewImage,2)/2 size(OverviewImage,2)/2],...
        'Ydata',[-size(OverviewImage,1)/2 size(OverviewImage,1)/2],...
        'FillValues',128);
    imwrite(OverviewImage_rotated_shifted,OverviewImageAlignedFileNameStr,'Compression','none');
    
  
end
matlabpool CLOSE;
for i=1:length(Files)
    OverviewAlignedDataFileNameStr =[secOverviewFiles.alignedDir{i} 'SectionOverviewAligned_' num2str(labels(i)) '.mat'];

    M=Ms{i};
    AlignmentParameters.r_offset = M(3,2);
    AlignmentParameters.c_offset = M(2,2);
    AlignmentParameters.AngleOffsetInDegrees = atan2(M(1,2),M(1,1))*180/pi;
    save(OverviewAlignedDataFileNameStr, 'AlignmentParameters'); 
end


disp('done');


