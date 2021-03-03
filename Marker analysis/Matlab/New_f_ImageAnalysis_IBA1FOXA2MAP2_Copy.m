%function  [ObjectsThisOrganoid] = FINAL_f_ImageAnalysisPerOperettaOrganoid_cell_count(Label, ch1, ch2, ch3, ch4, ChannelNames, PreviewPath);
function  [ObjectsThisOrganoid] = New_f_ImageAnalysis_IBA1FOXA2MAP2_Copy(Label, ch1, ch2, ch3, ch4, ChannelNames, PreviewPath);

%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    % vol(ch1, 0, 3000) % Alexa 488 >>> IBA1
    % vol(ch2, 0, 1000) % Alexa 647 >>> MAP2
    % vol(ch3, 0, 3000) % HOECHST 33342 >>> Hoechst imtool(max(ch2, [], 3))
    % vol(ch4, 0, 500) % TRITC >>> FOXA2

    %% Initialize variables
    NucleiMask = [];
    MAP2Mask = [];
    FOXA2Mask = [];
    IBA1Mask = []; 
    
    %% Segment nuclei
    %vol(ch3, 0, 2000)
    ch3BlurSmall = imfilter(double(ch3), fspecial('gaussian', 21, 1), 'symmetric');%vol(ch3BlurSmall)
    ch3BlurBig = imfilter(double(ch3), fspecial('gaussian', 21, 3), 'symmetric');%vol(ch3BlurBig) %%kind of flatfield corraction, to account for different bk in the pic
    ch3DoG = ch3BlurSmall - ch3BlurBig; %vol(ch3DoG, 0, 500, 'hot')
    NucleiMask = ch3DoG > 40; %vol(NucleiMask)
    NucleiMask = bwareaopen(NucleiMask, 200);%vol(NucleiMask)
    %bwareaopen excludes from the mask what is smaller than the specified
    %number of pixels
    ch3LP = imfilter(ch3, fspecial('gaussian', 11, 1), 'symmetric');%vol(ch3LP, 0, 10000, 'hot')
    NucMaskHigh =  (ch3LP > 2000) .* NucleiMask; %vol(NucMaskHigh, 0, 1) %%before3500 
    % NucMaskHigh discards what is dead
    NucMaskAlive = NucleiMask & ~NucMaskHigh; % vol(NucMaskAlive)
     
    %%  MAP2 (ch2)
  
    ch2MedFilt = []; 
    SizeZ = size(ch2, 3);
    parfor p = 1:SizeZ
        ch2MedFilt(:,:,p) = medfilt2(ch2(:,:,p));
    end
    %vol(ch2MedFilt, 0, 1000, 'hot')
     MAP2Mask = ch2MedFilt > 100; % vol(MAP2Mask)
    MAP2DoG = imfilter(ch2, fspecial('gaussian', 11, 1), 'symmetric') - imfilter(ch2, fspecial('gaussian', 31, 10), 'symmetric');
    %vol(MAP2DoG, 0, 20, 'hot')
    MAP2DoGMask = MAP2DoG > 10;
    %vol(MAP2DoGMask)
    MAP2Mask = MAP2Mask & MAP2DoGMask;
    %vol(MAP2Mask, 0, 1)
    %it(max(MAP2Mask, [], 3))
    %MAP2Mask = MAP2Mask & ~NucleiMask; %%%%%% IF THIS IS ON WILL CREATE AN ERROR IN THE %TUJ1
    MAP2Mask = bwareaopen(MAP2Mask, 200);
    %vol(MAP2Mask)
    

%% FOXA2 (ch4)
 
    ch4BlurSmall = imfilter(double(ch4), fspecial('gaussian', 21, 1), 'symmetric');%vol(ch4BlurSmall)
    ch4BlurBig = imfilter(double(ch4), fspecial('gaussian', 21, 3), 'symmetric');%vol(ch4BlurBig) %%kind of flatfield corraction, to account for different bk in the pic
    ch4DoG = ch4BlurSmall - ch4BlurBig; %vol(ch4DoG, 0, 50, 'hot')
    FOXA2Mask = ch4DoG > 15; %vol(FOXA2Mask)
    FOXA2Mask = bwareaopen(FOXA2Mask, 350);%vol(FOXA2Mask)
    %bwareaopen excludes from the mask what is smaller than the specified
    %number of pixels
    
  
   %%  IBA1 (ch1)
   
   ch1MedFilt = []; 
    SizeZ = size(ch1, 3);
    parfor p = 1:SizeZ
        ch1MedFilt(:,:,p) = medfilt2(ch1(:,:,p));
    end
    %vol(ch1MedFilt, 0, 3000, 'hot')
    IBA1Mask = ch1MedFilt > 200; % vol(IBA1Mask)
    IBA1DoG = imfilter(ch1, fspecial('gaussian', 75, 5) - fspecial('gaussian', 75, 25), 'symmetric');
    %IBA1DoG = gather(imfilter(gpuArray(ch1), gpuArray(fspecial('gaussian', 75, 5) - fspecial('gaussian', 75, 25)), 'symmetric'));
    %vol(IBA1DoG, 0, 50, 'hot')
    IBA1DoGMask = IBA1DoG > 20;
    %vol(IBA1DoGMask)
    
    IBATopHat = imfilter(imtophat(ch1, strel('disk', 15)), fspecial('gaussian', 11, 5), 'symmetric');%vol(IBATopHat, 0, 500)
    IBA1TopHatMask = IBATopHat > 150; % vol(IBA1TopHatMask)
    
    IBAStd = medfilt3(stdfilt(ch1, ones(3)));% vol(IBAStd, 0, 100)
    IBAStdMask = IBAStd > 55;
    IBAStdMask = bwareaopen(IBAStdMask,1000); %vol(IBAStdMask)
   
   
    IBA1Mask = imreconstruct(IBAStdMask,IBA1Mask & IBA1DoGMask & IBA1TopHatMask);
    %vol(IBA1Mask, 0, 1)
    %it(max(IBA1Mask, [], 3))
    %IBA1Mask = IBA1Mask & ~NucleiMask; %%%%%% IF THIS IS ON WILL CREATE ERROR IN THE %IBA1
    IBA1Mask = bwareaopen(IBA1Mask, 1000);
    %vol(IBA1Mask))
  
    %% Previews with mask
    
% Scalebar
    imSize = [size(ch1, 1), size(ch1, 2)];
    [BarMask, BarCenter] = f_barMask(200, 0.42, imSize, imSize(1)-200, 200, 25);
    %it(BarMask)
    
    PreviewHoechst = imoverlay2(imadjust(max(ch3,[],3),[0 0.075]), bwperim(max(NucleiMask,[],3)), [1 0 0]);
    PreviewHoechst = imoverlay2(PreviewHoechst, BarMask, [1 1 1]);
    % imtool(PreviewHoechst)
    
    PreviewNucMaskAlive = imoverlay2(imadjust(max(ch3,[],3),[0 0.07]), bwperim(max(NucMaskAlive,[],3)), [1 0 0]);
    PreviewNucMaskAlive = imoverlay2(PreviewNucMaskAlive, BarMask, [1 1 1]);
    % imtool(PreviewNucMaskAlive)

    PreviewIBA1 = imoverlay2(imadjust(max(ch1,[],3),[0 0.015]), bwperim(max(IBA1Mask,[],3)), [0 1 0]);
    PreviewIBA1 = imoverlay2(PreviewIBA1, BarMask, [1 1 1]);
    %imtool(PreviewIBA1)
    
    PreviewFOXA2 = imoverlay2(imadjust(max(ch4,[],3),[0 0.005]), bwperim(max(FOXA2Mask,[],3)), [1 0 0]);
    PreviewFOXA2 = imoverlay2(PreviewFOXA2, BarMask, [1 1 1]);
    %imtool(PreviewFOXA2)
    
    PreviewMAP2 = imoverlay2(imadjust(max(ch2, [], 3), [0 0.005]), bwperim(max(MAP2Mask,[],3)), [0 0 1]);
    PreviewMAP2 = imoverlay2(PreviewMAP2, BarMask, [1 1 1]);
    %imtool( PreviewMAP2)
    
    IdentityString = [Label.AreaName{:}, '_Idx_', num2str(Label.Idx)];
    imwrite(PreviewFOXA2, [PreviewPath, filesep, IdentityString, '_', 'FOXA2', '.png'])
    imwrite(PreviewHoechst, [PreviewPath, filesep, IdentityString, '_', 'Hoechst', '.png'])
    imwrite(PreviewIBA1, [PreviewPath, filesep, IdentityString, '_', 'IBA1', '.png'])
    imwrite(PreviewMAP2, [PreviewPath, filesep, IdentityString, '_', 'MAP2Mask', '.png'])
    imwrite(PreviewNucMaskAlive, [PreviewPath, filesep, IdentityString, '_', 'NucMaskAlive', '.png'])
       %% Feature extraction
    
    ObjectsThisOrganoid = table();
    ObjectsThisOrganoid.LabelIdx = {Label.Idx};
    ObjectsThisOrganoid.AreaName = {Label.AreaName};
    ObjectsThisOrganoid.NucMaskSum = sum(NucleiMask(:));
    ObjectsThisOrganoid.NucMaskAlive = sum(NucMaskAlive(:));
    ObjectsThisOrganoid.NucMaskHigh = sum(NucMaskHigh(:));
    ObjectsThisOrganoid.MAP2MaskSum = sum(MAP2Mask(:));
    ObjectsThisOrganoid.MAP2ByNucAlive = sum(MAP2Mask(:)) / sum(NucMaskAlive(:));
    ObjectsThisOrganoid.FOXA2MaskSum = sum(FOXA2Mask(:));
    ObjectsThisOrganoid.FOXA2ByNucAlive = sum(FOXA2Mask(:)) / sum(NucMaskAlive(:));
    ObjectsThisOrganoid.IBA1MaskSum = sum(IBA1Mask(:));
    ObjectsThisOrganoid.IBA1ByNucAlive = sum(IBA1Mask(:)) / sum(NucMaskAlive(:));
    ObjectsThisOrganoid.DeadCells = sum(NucMaskHigh(:))/sum(NucleiMask(:));
    ObjectsThisOrganoid.LiveNucCellCount = sum(NucMaskAlive(:)/67.9);
    ObjectsThisOrganoid.MAP2CellCount = sum(MAP2Mask(:)/404.3);
    ObjectsThisOrganoid.MAP2ByLiveCells = sum(ObjectsThisOrganoid.MAP2CellCount(:))/ (ObjectsThisOrganoid.LiveNucCellCount(:));
    ObjectsThisOrganoid.FOXA2CellCount = sum(FOXA2Mask(:)/67.9);
    ObjectsThisOrganoid.FOXA2ByLiveCells = sum(ObjectsThisOrganoid.FOXA2CellCount(:))/ (ObjectsThisOrganoid.LiveNucCellCount(:));
    ObjectsThisOrganoid.IBA1CellCount = sum(IBA1Mask(:)/146.6);
    ObjectsThisOrganoid.IBA1ByLiveCells = sum(ObjectsThisOrganoid.IBA1CellCount(:))/ (ObjectsThisOrganoid.LiveNucCellCount(:));
   
    

end

