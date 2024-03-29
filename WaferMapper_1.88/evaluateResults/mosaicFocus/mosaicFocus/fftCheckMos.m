%% Turn mosaics from waffer into viewable (centers and downsample) mosaics

%% Get Waffer folder information
wif = GetMyWafer;
%% View stitched

%% Parse XML
[tree, rootname, dom]=xml_read(wif.xml{end});
overlap = tree.MosaicSetup.TileOverlapXum;
xs = tree.MosaicSetup.TileWidth;
ys = tree.MosaicSetup.TileHeight;
res = tree.MosaicSetup.PixelSize;
pixelOverlap = fix(overlap / res * 1000);
minXfield = tree.MosaicSetup.Width;
minYfield = tree.MosaicSetup.Height;

%% Target Dir 
TPN = wif.dir; TPN = [TPN(1:end-1) 'Shaped\'];
TPNsav = [TPN 'quality\'];
if ~exist(TPNsav),mkdir(TPNsav);end

%% test Focus
colormap gray(256)
isize = 200;
subReg = [round(xs/2 - isize/2) round(xs/2 - isize/2) + isize-1];


for s = 4%1 : length(wif.sec) % run sections
    sprintf('reading tile %d of %d',s,length(wif.sec))
   
    %[tree, rootname, dom]=xml_read(wif.sec(s).xml);
    rc = wif.sec(s).rc;
    mosDim = max(rc,[],1);
    mos = zeros(mosDim*isize,'uint8');
    for t = 1:length(wif.sec(s).tile)
        
        ystart = (rc(t,1)-1)*isize+1;
        xstart = (rc(t,2)-1)*isize+1;
        I = imread(wif.sec(s).tile{t},'PixelRegion',{subReg,subReg});
        %[mfIs, listYX] =funSubFFT(I);

        subI = double(I);
        subI = subI - mean(subI(:));
        fftIy = abs(fft(subI,[],1));
        fftIx = abs(fft(subI,[],2));
        mfIy = mean(fftIy,2);
        mfIx = mean(fftIx,1);
        mfIy = mfIy(2:fix(length(mfIy)/2));
        mfIx = mfIx(2:fix(length(mfIx)/2));
        
        subplot(2,1,1)
        image(I)
        subplot(2,1,2)
        plot(mfIy),hold on
        plot(mfIx,'r'),hold off
        pause
        
        bgFFT = mean(mfI(end-10:end));
        mfI = mfI - bgFFT;
        mfIs(c,:) = mfI;


        



        
        
        
        mos(ystart:ystart+size(I,1)-1,xstart:xstart+size(I,2)-1) = 255-I;
        %image(mos),pause(.1)
    end
    image(mos),pause(1)
    % imwrite(mosB,[TPNsav wif.secNam{s} '.tif'],'Compression','none')
end










