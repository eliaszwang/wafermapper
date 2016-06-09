function [points,PixelRegion]=getSURFpointsfromFiles(Files,Options)
    %use OpenSurf to get SURF points from each section in the cell array of
    %full file paths
    %Options is a struct which contains the options for the OpenSurf points
    %see help OpenSurf for details
    %returns points, a 1xZ cell array containing the SURF points
   
    if ~exist('Options','var')
        Options.verbose=false;
        Options.init_sample=2;
        Options.octaves=3;
        Options.tresh=.001;
        Options.center_frac=1.0;
        Options.refsection=round(length(Files)/2);
    end
    if ~isfield(Options,'verbose')
        Options.verbose=0;
    end
    if ~isfield(Options,'init_sample')
         Options.init_sample=2;
    end
    if ~isfield(Options,'octaves')
        Options.octaves=3;
    end
    if ~isfield(Options,'tresh')
        Options.tresh=.001;
    end
    if ~isfield(Options,'center_frac')
        Options.center_frac=1.0;
    end
    if ~isfield(Options,'refsection')
        Options.refsection=round(length(Files)/2);
    end
    imagepath=Files{1};
    theinfo=imfinfo(imagepath);
    N=theinfo.Height;
    M=theinfo.Width;
    Z=length(Files);
    centerX=M/2;
    centerY=N/2;
    Cols=round([centerX-(Options.center_frac*M/2) centerX+(Options.center_frac*M/2)]);
    Rows=round([centerY-(Options.center_frac*N/2) centerY+(Options.center_frac*N/2)]);
    PixelRegion={Rows,Cols};
    readRegion=Options.center_frac<1.0; 
    %% calculate points of interest
    points=cell(1,Z);
    
    thesize=matlabpool('size');
    if thesize==0
        matlabpool OPEN;
    end
    SURFOptions.verbose=Options.verbose;
    SURFOptions.init_sample=Options.init_sample;
    SURFOptions.octaves=Options.octaves;
    SURFOptions.tresh=Options.tresh;
  
%     if (readRegion)
%         data=imread(Files{Options.refsection},'PixelRegion',PixelRegion);
%     else
%         data=imread(Files{Options.refsection}); 
%     end
%     data = PreFilterImage(data,Options);
%     refpoints=OpenSurf(data,SURFOptions);
%     numrefpoints=length(refpoints);
    
    parfor i=1:Z
        if (readRegion)
            data=imread(Files{i},'PixelRegion',PixelRegion);
        else
            data=imread(Files{i}); 
        end
        data = PreFilterImage(data,Options);
       
        points{i}=OpenSurf(data,SURFOptions);
            
%         if length(points{i})<.75*numrefpoints
%             disp('retaking');
%             disp([length(points{i}) numrefpoints]);
%             reSURFOptions=SURFOptions;
%             reSURFOptions.tresh=reSURFOptions.tresh*(max(length(points{i},1))/numrefpoints); %needs better threshold
%             points{i}=OpenSurf(data,reSURFOptions);
%             disp('after retake');
%             disp([length(points{i}) numrefpoints]);
%         end

        points{i}=shiftSURFpoints(points{i},-size(data,2)/2,-size(data,1)/2);

        disp([i Z]);
    end
    matlabpool CLOSE;
  
end
