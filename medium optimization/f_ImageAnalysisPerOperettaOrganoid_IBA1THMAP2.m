%function  [ObjectsThisOrganoid] = FINAL_f_ImageAnalysisPerOperettaOrganoid_cell_count(Label, ch1, ch2, ch3, ch4, ChannelNames, PreviewPath);
function  [ObjectsThisOrganoid] = New_f_ImageAnalysisPerOperettaOrganoid_IBA1THTUJ1(Label, ch1, ch2, ch3, ch4, ChannelNames, PreviewPath);

%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    % vol(ch1, 0, 3000) % Alexa 488 >>> TH
    % vol(ch2, 0, 1000) % Alexa 647 >>> MAP2
    % vol(ch3, 0, 3000) % HOECHST 33342 >>> Hoechst imtool(max(ch2, [], 3))
    % vol(ch4, 0, 500) % TRITC >>> IBA1

    %% Initialize variables
    NucleiMask = [];
    THMask = [];
    MAP2Mask = [];
    IBA1Mask = [];
    
    %% Segment nuclei
    %vol(ch3, 0, 2000)
    ch3BlurSmall = imfilter(double(ch3), fspecial('gaussian', 21, 1), 'symmetric');%vol(ch3BlurSmall)
    ch3BlurBig = imfilter(double(ch3), fspecial('gaussian', 21, 3), 'symmetric');%vol(ch3BlurBig) %%kind of flatfield corraction, to account for different bk in the pic
    ch3DoG = ch3BlurSmall - ch3BlurBig; %vol(ch3DoG, 0, 500, 'hot')
    NucleiMask = ch3DoG > 40; %vol(NucleiMask)
    NucleiMask = bwareaopen(NucleiMask, 500);%vol(NucleiMask)
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
     MAP2Mask = ch2MedFilt > 300; % vol(MAP2Mask)
    MAP2DoG = imfilter(ch2, fspecial('gaussian', 11, 1), 'symmetric') - imfilter(ch2, fspecial('gaussian', 31, 10), 'symmetric');
    %vol(MAP2DoG, 0, 20, 'hot')
    MAP2DoGMask = MAP2DoG > 40;
    %vol(MAP2DoGMask)
    MAP2Mask = MAP2Mask & MAP2DoGMask;
    %vol(MAP2Mask, 0, 1)
    %it(max(MAP2Mask, [], 3))
    %MAP2Mask = MAP2Mask & ~NucleiMask; %%%%%% IF THIS IS ON WILL CREATE AN ERROR IN THE %TUJ1
    MAP2Mask = bwareaopen(MAP2Mask, 200);
    %vol(MAP2Mask)   

% TH (ch1)
    ch1MedFilt = []; 
    SizeZ = size(ch1, 3);
    parfor p = 1:SizeZ
        ch1MedFilt(:,:,p) = medfilt2(ch1(:,:,p));
    end
    %vol(ch1MedFilt, 0, 1500, 'hot')
    THMask = ch1MedFilt > 50; % vol(THMask)
    THDoG = imfilter(ch1, fspecial('gaussian', 11, 1), 'symmetric') - imfilter(ch1, fspecial('gaussian', 31, 10), 'symmetric');
    %vol(THDoG, 0, 300, 'hot')
    THDoGMask = THDoG > 20;
    %vol(THDoGMask)
    THMask = THMask & THDoGMask;
    %vol(THMask, 0, 1)
    %it(max(THMask, [], 3))
    %THMask = THMask & ~NucleiMask; %%%%%% IF THIS IS ON WILL CREATE ERROR IN THE %TH  
    THMask = bwareaopen(THMask, 350);        
 
        %% TH Fragmentation
    
    % Define structuring element for surface detection
    Conn6 = strel('sphere', 1); % 6 connectivity
    % Detect surface
    SurfaceTH = THMask & ~(imerode(THMask, Conn6));
    %vol(SurfaceTH)
  
 %%  IBA1 (ch4)
   
    ch4MedFilt = []; 
    SizeZ = size(ch4, 3);
    parfor p = 1:SizeZ
        ch4MedFilt(:,:,p) = medfilt2(ch4(:,:,p));
    end
    %vol(ch41MedFilt, 0, 3000, 'hot')
    IBA1Mask = ch4MedFilt > 50; % vol(IBA1Mask)
    IBA1DoG = imfilter(ch4, fspecial('gaussian', 75, 5) - fspecial('gaussian', 75, 25), 'symmetric');
    %IBA1DoG = gather(imfilter(gpuArray(ch1), gpuArray(fspecial('gaussian', 75, 5) - fspecial('gaussian', 75, 25)), 'symmetric'));
    %vol(IBA1DoG, 0, 50, 'hot')
    IBA1DoGMask = IBA1DoG > 30;
    %vol(IBA1DoGMask)
    
    IBATopHat = imfilter(imtophat(ch4, strel('disk', 15)), fspecial('gaussian', 11, 5), 'symmetric');%vol(IBATopHat, 0, 500)
    IBA1TopHatMask = IBATopHat > 40; % vol(IBA1TopHatMask)
    
    IBAStd = medfilt3(stdfilt(ch4, ones(3)));% vol(IBAStd, 0, 100)
    IBAStdMask = IBAStd > 40;
    IBAStdMask = bwareaopen(IBAStdMask,1000); %vol(IBAStdMask)
  
    IBA1Mask = imreconstruct(IBAStdMask,IBA1Mask & IBA1DoGMask & IBA1TopHatMask);
    %vol(IBA1Mask, 0, 1)
    %it(max(IBA1Mask, [], 3))
    %IBA1Mask = IBA1Mask & ~NucleiMask; %%%%%% IF THIS IS ON WILL CREATE ERROR IN THE %IBA1
    IBA1Mask = bwareaopen(IBA1Mask,1000);
    %vol(IBA1Mask))


    %% Previews with mask
    
% Scalebar
    imSize = [size(ch1, 1), size(ch1, 2)];
    [BarMask, BarCenter] = f_barMask(200, 0.42, imSize, imSize(1)-200, 200, 25);
    %it(BarMask)
    
    PreviewHoechst = imoverlay2(imadjust(max(ch3,[],3),[0 0.09]), bwperim(max(NucleiMask,[],3)), [1 0 0]);
    PreviewHoechst = imoverlay2(PreviewHoechst, BarMask, [1 1 1]);
    % imtool(PreviewHoechst)
    
    PreviewNucMaskAlive = imoverlay2(imadjust(max(ch3,[],3),[0 0.09]), bwperim(max(NucMaskAlive,[],3)), [1 0 0]);
    PreviewNucMaskAlive = imoverlay2(PreviewNucMaskAlive, BarMask, [1 1 1]);
    % imtool(PreviewNucMaskAlive)

    PreviewIBA1 = imoverlay2(imadjust(max(ch4,[],3),[0 0.01]), bwperim(max(IBA1Mask,[],3)), [0 1 0]);
    PreviewIBA1 = imoverlay2(PreviewIBA1, BarMask, [1 1 1]);
    %imtool(PreviewIBA1)
    
    PreviewTH = imoverlay2(imadjust(max(ch1,[],3),[0 0.01]), bwperim(max(THMask,[],3)), [1 0 0]);
    PreviewTH = imoverlay2(PreviewTH, BarMask, [1 1 1]);
    %imtool(PreviewTH)
    
    PreviewMAP2 = imoverlay2(imadjust(max(ch2, [], 3), [0 0.02]), bwperim(max(MAP2Mask,[],3)), [0 0 1]);
    PreviewMAP2 = imoverlay2(PreviewMAP2, BarMask, [1 1 1]);
    %imtool( PreviewMAP2)
    
    IdentityString = [Label.AreaName{:}, '_Idx_', num2str(Label.Idx)];
    imwrite(PreviewTH, [PreviewPath, filesep, IdentityString, '_', 'TH', '.png'])
    imwrite(PreviewHoechst, [PreviewPath, filesep, IdentityString, '_', 'Hoechst', '.png'])
    imwrite(PreviewIBA1, [PreviewPath, filesep, IdentityString, '_', 'IBA1', '.png'])
    imwrite(PreviewMAP2, [PreviewPath, filesep, IdentityString, '_', 'MAP2Mask', '.png'])
    imwrite(PreviewNucMaskAlive, [PreviewPath, filesep, IdentityString, '_', 'NucMaskAlive', '.png'])

chEmpty = zeros(size(ch1),'uint16');
THRed = cat(3,imadjust(max(ch1, [], 3)), max(chEmpty,[],3), max(chEmpty,[],3)); %imtool(THRed)
IBAGreen = cat(3, max(chEmpty,[],3), imadjust(max(ch4, [], 3)),max(chEmpty,[],3)); %imtool(IBAGreen)
NucBlue = cat(3, max(chEmpty,[],3), max(chEmpty,[],3),imadjust(max(ch3, [], 3))); %imtool(NucBlue)
MAP2white = imoverlay2(imadjust(max(ch1,[],3),[0 0.02]), [1 1 1]); %imtool(MAP2white)
    %% Feature extraction
    
    ObjectsThisOrganoid = table();
    ObjectsThisOrganoid.LabelIdx = {Label.Idx};
    ObjectsThisOrganoid.AreaName = {Label.AreaName};
    ObjectsThisOrganoid.NucMaskSum = sum(NucleiMask(:));
    ObjectsThisOrganoid.NucMaskAlive = sum(NucMaskAlive(:));
    ObjectsThisOrganoid.DeadCells = sum(NucMaskHigh(:))/sum(NucleiMask(:));
    ObjectsThisOrganoid.NucMaskHigh = sum(NucMaskHigh(:));
    ObjectsThisOrganoid.MAP2MaskSum = sum(MAP2Mask(:));
    ObjectsThisOrganoid.MAP2ByNucAlive = sum(MAP2Mask(:)) / sum(NucMaskAlive(:));
    ObjectsThisOrganoid.THMaskSum = sum(THMask(:));
    ObjectsThisOrganoid.THByMAP2 = sum(THMask(:)) / sum(MAP2Mask(:));
    ObjectsThisOrganoid.THByNucAlive = sum(THMask(:)) / sum(NucMaskAlive(:));
    ObjectsThisOrganoid.THFragmentation = sum(SurfaceTH(:)) / sum(THMask(:));
    ObjectsThisOrganoid.THByMAP2 = sum(THMask(:)) / sum(MAP2Mask(:));
    ObjectsThisOrganoid.THByNucAlive = sum(THMask(:)) / sum(NucMaskAlive(:));
    ObjectsThisOrganoid.IBA1MaskSum = sum(IBA1Mask(:));
    ObjectsThisOrganoid.IBA1ByNucAlive = sum(IBA1Mask(:)) / sum(NucMaskAlive(:));

end


    
   
  
