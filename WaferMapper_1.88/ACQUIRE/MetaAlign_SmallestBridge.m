function [finalTransforms,num_inliers,notSectionFlag]=MetaAlign_SmallestBridge(points,Options,Files,MatFiles)
%Takes a set of point features across Z sections (currently written to use OPENSurf
%(points is thus a 1xZ cell array)
%taken from the virtual tiff stack located in StackDirectory (currently i read in
%the images so one could display the matches, though it isn't strictly
%necessary for the calculation
%
%from these set of features, it calculates a set of similarity
%tranformations that
if ~exist('Options','var')
    Options.verbose=0;
    Options.Nbest=150;
    Options.dist_thresh=.02; 
    Options.ref_sect=1;
    Options.det_thresh=.25;
    Options.min_inliers=20;
    Options.max_dist=15;
end
if ~isfield(Options,'verbose')
    Options.verbose=0;
end
if ~isfield(Options,'Nbest')
    Options.Nbest=150;
    % This specifies that it should find at most the 50 best points of correspondance
    % between any two sections before running ransac to find the
    % set of points amongst those which are best explained by a
    % single rigid transformation
end
if ~isfield(Options,'dist_thresh')
    Options.dist_thresh=.05; %this is how close in percentage of the size of the image read
    %the transformed points have to meet their matched points in order
    %to be counted as an inlier
end
if ~isfield(Options,'ref_sect')
    Options.ref_sect=1; %This is how much the determinant can be different than 1 and still be a valid transform for ransac
end
if ~isfield(Options,'det_thresh')
    Options.det_thresh=.25; %This is how much the determinant can be different than 1 and still be a valid transform for ransac
end
if ~isfield(Options,'min_inliers')
     Options.min_inliers=20; %This is how many inliers are necessary for it to be considered a good link
end
if ~isfield(Options,'max_dist')
     Options.max_dist=15; %This is how many inliers are necessary for it to be considered a good link
end

theinfo=imfinfo(Files{1}); %pull out the info from the first section
N=theinfo.Height; %this is the Height/number of rows
M=theinfo.Width; %this is the Width/number of columns
Z=length(Files); %this is the number of sections in stack
%useful for debugging

%setup the ransac options
ransacCoef.minPtNum=2;   %for a rigid or similarity transform 2 is the number needed
ransacCoef.iterNum= Options.Nbest; %run through a number of iterations equal
% to the number of points
if isfield(Options,'PixelRegion')
    Cols=Options.PixelRegion{1};
    ransacCoef.thDist= Options.dist_thresh*(Cols(2)-Cols(1));
else
    ransacCoef.thDist= Options.dist_thresh*M;
end
ransacCoef.thInlrRatio=.01; %at least 1 percent should be right
ransacCoef.thDet=Options.det_thresh;

%setup tracking of shortest good Link
transformDepth=-ones(1,Z);
transformedTo=zeros(1,Z);
transforms=cell(1,Z);
Pos1s=cell(1,Z);
Pos2s=cell(1,Z);
Inliers=cell(1,Z);
num_inliers=zeros(1,Z);
notSectionFlag=zeros(1,Z);

transformDepth(Options.ref_sect)=0;

[dest,source]=meshgrid(1:Z,1:Z);

dist_source_to_dest=abs(source-dest);
level=0;
changed=1;
Pos1s_ij=cell(Z,Z);
Pos2s_ij=cell(Z,Z);
inliers_ij=cell(Z,Z);
transforms_ij=cell(Z,Z);
numinliers_ij=zeros(Z,Z);

while and(sum(transformDepth==-1)>0,changed==1) %continue until all the sections have been linked, or nothing was done
    changed=0; %reset the flag to check that this level did anything
    ToConnect=find(transformDepth==-1); %list out who still needs to be linked up
    goodDests=find(transformDepth==level); %list out the list of sections one could link to at this level
    
    for source_ind=1:length(ToConnect) %loop over all the sections we need to link
        source=ToConnect(source_ind); %pull out the section we are linking
        distances=dist_source_to_dest(source,goodDests); %pull out the distances from this section to all the valid destinations
        [sorteddistances,order]=sort(distances,2,'ascend'); %sort them in ascending order
        close_enough=find(sorteddistances<=Options.max_dist); %get the list of indices into this order of distances that are close enough
        destinations=goodDests(order(close_enough)); %get the list of destinations in order
        for dest_ind=1:length(destinations) %loop over the destinations
            dest=destinations(dest_ind); %pull out the destination
            if and(exist('MatFiles','var'),isfield(Options,'PixelRegion')) %remove the overlapping regions of the images if the info was passed
                [points1,points2]=RemovePointsFromOverlapZone(points{source},points{dest},Options.PixelRegion,MatFiles{source},MatFiles{dest});
            else %otherwise just pull out the default set of points
                points1=points{source};
                points2=points{dest};
            end
           
            % return the Options.Nbest SURF matches between these two sections
            [Pos1,Pos2]=find_best_SURF_match(points1,points2,Options.Nbest);
            
            %now use those correspondances using a rigid transformation
            %need to transpose list of points to make them  (Y,X)x Nbest, and then flipdim to make it (X/Y)xNbest
            
            %run Ransac
            [f inlierIdx] = ransac1( flipdim(Pos1',1),flipdim(Pos2',1),ransacCoef,@fit_rigid,@EuclideanDistance,@getDeterminantOfAffineTransform);
            
            %keep track of this fit, just in case we need it later
            Pos1s_ij{source,dest}=Pos1;
            Pos2s_ij{source,dest}=Pos2;
            inliers_ij{source,dest}=inlierIdx;
            transforms_ij{source,dest}=f;  
            numinliers_ij(source,dest)=length(inlierIdx);
            
            %if there are more than the minimum needed inliers than we will consider this a
            %match and record it as such
            if(length(inlierIdx)>=Options.min_inliers)
                 disp([source dest]); %display some progress
                changed=1; %flag this level as having accomplished something
                transformDepth(source)=level+1; %mark this link as being done at the appropriate level
                transformedTo(source)=dest; %note which section this was linked to
                transforms{source}=f; %save the transform
              
                Pos1s{source}=Pos1; %record the matching positions/tranforms
                Pos2s{source}=Pos2; 
           
                Inliers{source}=inlierIdx; %record the final inliers
                num_inliers(source)=length(inlierIdx); %record how many inliers in this link
                if Options.verbose                    
                    %to visualize the results
                    I1=imread(Files{source});
                    I2=imread(Files{dest});
                    figure(4);
                    clf;
                    hold off;
                    PlotCorrespondances(I1,I2,Pos1,Pos2,inlierIdx)        
                    figure(5);
                    clf;
                    I1t=imtransform(I1,transforms{source},'bicubic',...
                        'UData',[-size(I2,2)/2 size(I2,2)/2],...
                        'VData',[-size(I2,1)/2 size(I2,1)/2],...
                        'XData',[-size(I2,2)/2 size(I2,2)/2],...
                        'Ydata',[-size(I2,1)/2 size(I2,1)/2],...
                        'FillValues',128);
                    
                    I = zeros([size(I1,1) size(I1,2) 3],'double');
                    I(:,:,1)=double(I1t) / 255.0;
                    I(:,:,2)=double(I2)/ 255.0;
                    %display the resulting overlap in red/green
                
                    imshow(I); 
                    pause;
                    
                end
                %now we don't want to consider any more mathces for this
                %source, so break out of the loop.
                break;
            end
        end
    end
    level=level+1; %increment the level
end

%now we need to iterate through the sections
%to produce this final set of composed transforms
finalTransforms=cell(1,Z);


for i=1:Z
   
    %if this section got linked to another section
    if transformedTo(i)~=0
        %this is what should happen 99% of the time
        finalTransforms{i}=CreateComposedTranformFromLinkedList(i,transformedTo,transforms);
    else
        %if it isn't linked to anyone it is either
        if i==Options.ref_sect %the reference section
            %in which case we  should give it the identity transform
            finalTransforms{i}=maketform('affine',[1 0 0; 0 1 0; 0 0 1]);%then it should default to Identity transform
        else %or it's an orphan section
            %and we need to decide what to do with it
            numinlier_vec=numinliers_ij(i,:);
            
            %look through the fits that we did, and order them by number of
            %inliers
            [sorted_numinliers,sorted_order]=sort(numinlier_vec,2,'descend');
            
            %loop through these in sorted order, and ask user if they want
            %to accept one of these fits even with fewer inliers than
            %MinInliers.. or perhaps it isn't a section and should be noted
            %as such.
            goodones=find(sorted_numinliers>0);
            for j=1:length(goodones)
                %this is the candidate destination
                dest=sorted_order(goodones(j));
                source=i;
                %read in the images
                I1=imread(Files{source});
                I2=imread(Files{dest});
                figure(3);
                clf;
                hold off;
                PlotCorrespondances(I1,I2,Pos1s_ij{source,dest},Pos2s_ij{source,dest},inliers_ij{source,dest});
                %ask the user what they want to do
                choice=questdlg(['Orphan Section ' num2str(i)],'Orphan Section Dialog','Accept','See Another','Remove Section','Accept');
                close(3);
                switch choice
                    case 'Accept'
                        %if they accept this fit, update the linked list
                        %and make a final transform as if it was proper all
                        %along
                        transformedTo(i)=dest;
                        transforms{i}=transforms_ij{i,dest};
                        transformDepth(i)=transformDepth(dest)+1;
                        finalTransforms{i}=CreateComposedTranformFromLinkedList(i,transformedTo,transforms);
                        break; %and then break out of the loop
                    case 'See Another'
                        continue; %let the loop continue
                    case 'Remove Section'
                        notSectionFlag(i)=1; %note this as not a section
                        finalTransforms{i}=maketform('affine',[1 0 0; 0 1 0; 0 0 1]);%default to Identity transform
                        break;
                end
                
            end
            %we got through the loop or we broke out
            if and(notSectionFlag==0,transformedTo(i)==0) %if we did not break out
               %then this section is still orphaned and will need to be
               %fixed manually, so alert the user
               errordlg(['Section ' num2str(i) ' file:' Files{i} ' Still Orphaned Must Fix Manually']);
               finalTransforms{i}=maketform('affine',[1 0 0; 0 1 0; 0 0 1]);%will default to Identity transform
            end
            
        end
    end
end


function finalTransform=CreateComposedTranformFromLinkedList(startsection,transformedTo,transforms)
if transformedTo(startsection)~=0
    TransformList=[];
    currInd=startsection;
    %follow the chain to the reference section
    while transformedTo(currInd)~=0
        TransformList=[transforms{currInd} TransformList]; %append this transform
        currInd=transformedTo(currInd); %update currInd to the next section in chain to reference
    end
    
    %make a composed transform if necessary
    if length(TransformList)>1
        finalTransform=maketform('composite',TransformList);
    elseif length(TransformList)==1 %otherwise it is just this one
        finalTransform=TransformList(1);
    end
else
    finalTransform=[];
end



