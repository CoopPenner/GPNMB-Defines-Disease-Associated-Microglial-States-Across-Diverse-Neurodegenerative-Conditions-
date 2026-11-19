

%% Chemotaxis figure maker

%WTVid='/Volumes/PC60/iMicroVids/setsPassedQC/PBS/WTKOScratchAssayPBSOnly_6_21_26_allChannels/WT/WT_Scratch_PBS_Well4 (1).mp4';

%vidAtPlay='/Volumes/PC60/iMicroVids/setsPassedQC/PBS/WTKOScratchAssayPBSOnly_6_21_26_allChannels/KO/KO_PBS_Scratch_Well4 (2).mp4';
%vidAtPlay='/Volumes/PC60/iMicroVids/setsPassedQC/PBS/WTKOScratchAssayPBSOnly_6_21_26_allChannels/KO/KO_PBS_Scratch_Well4 (1).mp4';
vidAtPlay='/Volumes/PC60/iMicroVids/setsPassedQC/PBS/WTKOScratchAssayPBSOnly_6_21_26_allChannels/WT/WT_Scratch_PBS_Well4 (2).mp4'

         vidObj = VideoReader(vidAtPlay);



             frameIdx = 1;

     numFrames = floor(vidObj.Duration * vidObj.FrameRate)-1;
        dapiFrames = cell(numFrames, 1);
        BodipyFrames= cell(numFrames,1);
        refFrames=cell(numFrames, 1);
        while hasFrame(vidObj)
            rgbFrame = readFrame(vidObj);
        
            % Extract the blue channel (DAPI)
            blueChannel = rgbFrame(:,:,3);  % Blue is the 3rd channel in RGB
            greenChannel= rgbFrame(:,:,2);
            dapiFrames{frameIdx} = blueChannel;
            BodipyFrames{frameIdx}=greenChannel;
            refFrames{frameIdx}=rgbFrame;
            frameIdx = frameIdx + 1;
        end
        

        % Grab first and last frames from the already-extracted cell arrays
firstFrameRGB = refFrames{1};
lastFrameRGB  = refFrames{92};


     totfirstIm = cat(3, nan(size(BodipyFrames{1})), BodipyFrames{1}, dapiFrames{1});     % 3D array HxWx3


% Display side by side
figure;
imshow(totfirstIm);
title('First Frame');

exportgraphics(gca, '/Volumes/PC60/iMicroVids/FinalFigures/jointAnalyses/PBS/WTFirstFrameExemplar2NoChangeFluoro.pdf', 'ContentType', 'vector')


     totLastIm = cat(3, nan(size(BodipyFrames{1})), BodipyFrames{end}, dapiFrames{end});     % 3D array HxWx3




     
figure
imshow(totLastIm);
title('Last Frame');

exportgraphics(gca, '/Volumes/PC60/iMicroVids/FinalFigures/jointAnalyses/PBS/WTLastFrameExemplar2NoChangeFluoro.pdf', 'ContentType', 'vector')
