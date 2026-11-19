function [cleanedMask,CC,nucFeatures] = lipidMaskMaker(frameAtPlay,minSize, maxSize, circCut,strelRadiusSize,strelType,backGroundCorrect,backGroundGaussian, demoPlot)



%convert to grayscale if necessary
        if size(frameAtPlay,3)==3
        frameAtPlay=rgb2gray(frameAtPlay);
        warning('running algorithm on rgb image, consider using only nuclear channel')
        end


if backGroundCorrect

I = im2double((frameAtPlay));
I = imflatfield(I, backGroundGaussian);   % radius in pixels; try 40–120
frame2use=I;

else
frame2use=frameAtPlay;
end

% Binarize the  first frame of the image using Otsu's method
% Binarize the  first frame of the image using Otsu's method
%nuclearMask = imbinarize(frame2use, 'global'); 
nuclearMask = imbinarize(frame2use, 'adaptive','Sensitivity',.5); 





% 
% figure
% imshowpair(frameAtPlay,nuclearMask,'montage')

se = strel(strelType,strelRadiusSize);  % now using strel to fill in small gaps between different parts of a given nuclei. 


nuclearMaskOpened = imopen(nuclearMask, se);  % remove specks first
nuclearMaskClosed = imclose(nuclearMaskOpened, se);  % then close gaps



% Apply initial size filtering
cleanedMask = bwareaopen(nuclearMaskClosed, minSize);

% figure
% imshowpair(frameAtPlay,cleanedMask,'montage')


% Get connected components and nuclear features
CC = bwconncomp(cleanedMask);
nucFeatures = regionprops("table", CC, ...
    "Area", "Centroid", "EquivDiameter", "Solidity", "Circularity");

% Extract features
areaAtPlay = nucFeatures.Area;
circ = nucFeatures.Circularity;

% Define  nuclei more likely to be artifactual
badIdx = areaAtPlay > maxSize | circ < circCut;
% how I calculate distance


% Set those components to false in the binary mask
labeledMask = labelmatrix(CC);
for ss = find(badIdx)'
    cleanedMask(labeledMask == ss) = 0;
end


% regenerate connected components


CC = bwconncomp(cleanedMask);

nucFeatures = regionprops("table", CC, ...
    "Area", "Centroid", "EquivDiameter", "Solidity", "Circularity","PixelIdxList","PixelList");








if demoPlot

figure

subplot(2,2,1)
imshow(frameAtPlay)

subplot(2,2,2)
imshow(nuclearMask)

subplot(2,2,3)
imshow(nuclearMaskClosed)

subplot(2,2,4)
imshow(cleanedMask)
end



end


%% initial cell selection... this will be the mask that we use to track

% my initial approach will be to prioritize eliminating low intensity
% nuclei and differentiating nearby cells. Ie removing false positives 

% avgFrame = (double(rgb2gray(frames{1})) + double(rgb2gray(frames{2})))  / 2;