%% DQ-BSA signal analysis



%% initializing relevant variables and loading video
addpath('/Users/pennerc/Documents/Documents - BGS-STU23-138/MATLAB/GPNMB_Global_Analysis/GPNMB-Global-Analysis');
addpath('/Users/pennerc/Documents/Documents - BGS-STU23-138/GitHub/GPNMB_GlobalANalysis_Updated')

% first upload relevant video

homePath= '/Volumes/PC60/iMicroVids/lysoFunction/7_17_26/KO_DQBSA';


genoTypePrefix='KO'; % I do this just as a safe check
imageTypePrefix='DQBSA';

outputSuffix='7_17_26';

frameRate=15; %one image every 15 minutes



%Mask Making factors


%nucleus (currently this isn't used since the hoescht was inconsistent) 


%this was the biggest issue for this pipeline, ultimately I got it working
%decently but opted to just go based on the BSA. My initial run through
%showed that BSA objects did not actually increase over time. So all cells
%that will become positive start positive. this allowed us to implement a
%relativley straight forward approach of simply looking at avg fluoro per
%image over time. SIngle cell tracking will be helpful to pick out
%individual dynamics but will be most beneifical with another marker to
%identify, for example, concurrent lipid flux or some such thing. In this
%case the BODIPY became totally blown out for the KO making that channel
%not usable. so it goes


minSize=200; % 200 for 10x 250 for 10x lmao this is the minimum size an object can be (some dapi debris always present)
maxSize=2500;   %original run 600 % sometimes multiple cells become clumped (as in a swarm) I just remove these
circCut=0.33; % originally did .05 this is a cutoff for the minimum circularity of a detected object, less circular objects tend to be debris n whatnot.  
detectSensitivity=.5; % for now this is meaningless insofar as we are using OTSU's this is the sensitivity of our initial binarization step. I've found .6 is a good middle ground but decrease if there are a lot of artifactual objects
strelRadiusSize=3; %this is the pixel radius of our structuring element. We first apply an erosion of all pixels to remove noise and then proceed with dilation procedure
%a larger strel Radius Size will cause objects that are further apart ( 4
%vs 5 vs 6 pixels) to become aggregated
demoPlot=false; %this will just create a figure with the masks at each of our cleanup stages
totPlot=false; %this will output a montage of the raw image vs the mask for each frame
dapiFluoroIncrease=10; %this was necessary to get good detection all fluoro in a given channel multiplied by this factor 



  %Bodipy Staining
    
    BSAMAxSize=30000; % 
    BSAMinSize=400; %  edges of processes sometimes separated by strel
    BSADetectSensitivity=.5;
    BSARadiusSize=5; %here strel type is always arbitrary
    BSACircCut=0;
    BSAStrelType='arbitrary';



% % trajectory filtering features

%these are hold overs for if we go single cell
% maxJump=30; %this is the upper limit of how much a nucleus can move in one frame
% confThreshold=.05; %this is the lower limit of our confidence in tracking a nucleus from one point to another.
% areaThreshold=1; %this is how much of a difference in size from frame to frame you'll tolerate this will be protection against fusion or splitting events
% 

%extract files

allFilesAtPlay = dir(fullfile(homePath, '*.mp4'));
allFilesAtPlay = allFilesAtPlay(~startsWith({allFilesAtPlay.name}, '._')); % mac keeps these weird ghost files with ._ prefix

allFilesAtPlay= allFilesAtPlay(  contains( {allFilesAtPlay.name}, genoTypePrefix)  & ...
(contains( {allFilesAtPlay.name}, imageTypePrefix)));




nucMaskCollect=cell(1,length(allFilesAtPlay)); nucFeatureCollect= cell(1,length(allFilesAtPlay));rawNucCollect= cell(1,length(allFilesAtPlay));
BSAMaskCollect=cell(1,length(allFilesAtPlay)); BSAFeatureCollect= cell(1,length(allFilesAtPlay)); BSARawCollect= cell(1,length(allFilesAtPlay));

allResults=table(); %this is where the per cell per frame fluoresence gets accumulated



%% start of for Loop through all vids

 for vidNum=1: length(allFilesAtPlay)            
        
        vidName=[homePath,'/', allFilesAtPlay(vidNum).name];
 
         %KO_Scratch_aB (1).mp4
         vidObj = VideoReader(vidName);
   
    
        % loading objects into frames
        
        disp('loading video data, this will take a few seconds')
        
        %preallocating objects
        
        numFrames = floor(vidObj.Duration * vidObj.FrameRate)-1;
        dapiFrames = cell(numFrames, 1);
        BSAFrames= cell(numFrames,1);
        refFrames=cell(numFrames, 1);
        frameIdx = 1;

        while hasFrame(vidObj)
            rgbFrame = readFrame(vidObj);
        
            % Extract the blue channel (DAPI)
            blueChannel = rgbFrame(:,:,3);  % Blue is the 3rd channel in RGB
            redChannel= rgbFrame(:,:,1);
            dapiFrames{frameIdx} = blueChannel;
            BSAFrames{frameIdx}=redChannel;
            refFrames{frameIdx}=rgbFrame;
            frameIdx = frameIdx + 1;
        end
        
%I'm not sure why but the last frame is always blanked out, have seen this
%before...


        
        dapiMaskCollect=cell(numFrames, 1); %we will collect both masks and info about them
        frameDataCollect=cell(numFrames, 1);
        objNumCollect=nan(1,numFrames);

        BSAMaskFrames=cell(numFrames, 1); %per frame BSA masks and features for this video
        BSAFeatureFrames=cell(numFrames, 1);
        
        
        
            if ~demoPlot
                    frames2Assess=length(BSAFrames);
                    disp('video loaded, generating masks for each frame of movie')
        
            elseif demoPlot
                frames2Assess=6;
                disp('outputting image of binarized masks')
                figure
        
            end
        
        
                for bb=1: frames2Assess
                
                       
                   [dapiMaskCollect{bb},CC,frameDataCollect{bb}] = pulseChaseMaskMakerDapi(dapiFrames{bb}*dapiFluoroIncrease,minSize, maxSize, circCut, detectSensitivity,strelRadiusSize,demoPlot);
                
                        
                        objNumCollect(bb)=height(frameDataCollect{bb});


      [BSAMaskFrames{bb},~,BSAFeatureFrames{bb}] = aBNuclearMaskMakerTrem2Specific(BSAFrames{bb},BSAMinSize, BSAMAxSize, BSACircCut,strelRadiusSize,BSAStrelType, BSADetectSensitivity, demoPlot);

                
%bodipFrames{bb}

     totIm = cat(3,  BSAFrames{bb}, nan(size(BSAFrames{bb})), dapiFrames{bb}*20);     % 3D array HxWx3
     maskImTot = cat(3,  BSAMaskFrames{bb}, nan(size(BSAFrames{bb})), dapiMaskCollect{bb}*100  );  
    if totPlot  
     figure
     imshowpair(totIm,maskImTot,'montage')
     title(["over view of  ", allFilesAtPlay(vidNum).name])
    end
                


% here let's just do a simple fluoresence extraction, for this frame I will
% do this as an initial test. If the number of objects don't decrease it
% should be fine... 

        BSAFeaturesAtPlay=BSAFeatureFrames{bb};
        rawBSAIm=double(BSAFrames{bb});
        bodipMaskAtPlay=BSAMaskFrames{bb};

            if ~isempty(BSAFeaturesAtPlay)

                %bodipFluoroCollectNuc=nan(1,height(nucFeaturesAtPlay));
                %%signal associated with the nucleus
                BSAFluoroCollectCell=nan(1,height(BSAFeaturesAtPlay)); %signal in the cell body
                BSAAreaCollectCell=nan(1,height(BSAFeaturesAtPlay)); %bodipy positive pixels in that ROI



                    BSAFluoroCollect=nan(1,height(BSAFeaturesAtPlay));
                    BSAFluoroStd=nan(1,height(BSAFeaturesAtPlay));
                    BSAFluoroMax=nan(1,height(BSAFeaturesAtPlay));
                    BSAFluoroMin=nan(1,height(BSAFeaturesAtPlay));
                    BSAFluoroInt=nan(1,height(BSAFeaturesAtPlay)); %integrated density, area*mean
                    
                    for nn=1:height(BSAFeaturesAtPlay)
                        pixVals=rawBSAIm(BSAFeaturesAtPlay.PixelIdxList{nn}); %literally take pixel values in detected region
                        BSAFluoroCollect(nn)=nanmean(pixVals); %extract mean
                        BSAFluoroStd(nn)=nanstd(pixVals);
                        BSAFluoroMax(nn)=max(pixVals);
                        BSAFluoroMin(nn)=min(pixVals);
                        BSAFluoroInt(nn)=nansum(pixVals);
                    end
                    
                    nObj=height(BSAFeaturesAtPlay);




                nCells=height(BSAFeaturesAtPlay);

        tableAtPlay=table( repmat(string(allFilesAtPlay(vidNum).name), nObj,1), ...
        repmat(string(genoTypePrefix), nObj,1), ...
        repmat(bb, nObj,1), ...
        repmat(bb*frameRate, nObj,1), ...
        (1:nObj)', ...
        BSAFeaturesAtPlay.Centroid(:,1), ...
        BSAFeaturesAtPlay.Centroid(:,2), ...
        BSAFeaturesAtPlay.Area, ...
        BSAFluoroCollect', ...
        BSAFluoroStd', ...
        BSAFluoroMin', ...
        BSAFluoroMax', ...
        BSAFluoroInt', ...
        'VariableNames', {'video','genotype','frame','timeMin','objInFrame','centroidX','centroidY','BSAArea','BSAFluoroMean','BSAFluoroStd','BSAFluoroMin','BSAFluoroMax','BSAFluoroInt'});
              
        
        
        allResults=[allResults; tableAtPlay];
%nice, easy peezy, but we want something better... let's shoot for single
%cell tracking. 


            end

                   
                    

                end


        %stash this video's per frame output in the top level collectors
        nucMaskCollect{vidNum}=dapiMaskCollect; nucFeatureCollect{vidNum}=frameDataCollect; rawNucCollect{vidNum}=dapiFrames;
        BSAMaskCollect{vidNum}=BSAMaskFrames; BSAFeatureCollect{vidNum}=BSAFeatureFrames; BSARawCollect{vidNum}=BSAFrames;

        fprintf('Done: %s (%d frames)\n', allFilesAtPlay(vidNum).name, frames2Assess);

 end


writetable(allResults, fullfile(homePath, sprintf('%s_%s_bodipyPerFrame_%s.csv', genoTypePrefix, imageTypePrefix, outputSuffix)));





%% summarize, one row per video per timepoint above we output on a per object basis... 
vidList=unique(allResults.video);
meanTable=table();
for vv=1:length(vidList)
    thisVid=allResults(allResults.video==vidList(vv), :);
    frameList=unique(thisVid.frame);
for ff=1:length(frameList)
        thisFrame=thisVid(thisVid.frame==frameList(ff), :);
        nObj=height(thisFrame);
        meanTable=[meanTable; table( vidList(vv), ...
            frameList(ff), ...
            unique(thisFrame.timeMin), ...
            nObj, ...
            nanmean(thisFrame.BSAFluoroMean), ...
            nanstd(thisFrame.BSAFluoroMean), ...
            nanstd(thisFrame.BSAFluoroMean)/sqrt(nObj), ...
            nanmean(thisFrame.BSAArea), ...
            nanstd(thisFrame.BSAArea)/sqrt(nObj), ...
            nansum(thisFrame.BSAArea), ...
            nanmean(thisFrame.BSAFluoroInt), ...
            nansum(thisFrame.BSAFluoroInt), ...
            nanmean(thisFrame.BSAFluoroMax), ...
'VariableNames', {'video','frame','timeMin','nObj','meanFluoro','stdFluoro','semFluoro','meanArea','semArea','totArea','meanIntDen','totIntDen','meanFluoroMax'})];
end
end
%% mean fluoresence per object, one line per video, SEM band
figure; hold on
cmap=lines(length(vidList));
for vv=1:length(vidList)
    thisVid=meanTable(meanTable.video==vidList(vv), :);
    xFill=[thisVid.timeMin; flipud(thisVid.timeMin)];
    yFill=[thisVid.meanFluoro+thisVid.semFluoro; flipud(thisVid.meanFluoro-thisVid.semFluoro)];
    fill(xFill, yFill, cmap(vv,:), 'FaceAlpha', .15, 'EdgeColor','none', 'HandleVisibility','off');
    plot(thisVid.timeMin, thisVid.meanFluoro, '-o', 'Color', cmap(vv,:), 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', vidList(vv));
end
xlabel('time (min)'); ylabel('mean DQ-BSA fluoresence per object')
title([genoTypePrefix, ' DQ-BSA intensity over time'])
legend('Location','eastoutside','Interpreter','none'); box off
%% mean object area, SEM band this didn't really work
figure; hold on
for vv=1:length(vidList)
    thisVid=meanTable(meanTable.video==vidList(vv), :);
    xFill=[thisVid.timeMin; flipud(thisVid.timeMin)];
    yFill=[thisVid.meanArea+thisVid.semArea; flipud(thisVid.meanArea-thisVid.semArea)];
    fill(xFill, yFill, cmap(vv,:), 'FaceAlpha', .15, 'EdgeColor','none', 'HandleVisibility','off');
    plot(thisVid.timeMin, thisVid.meanArea, '-o', 'Color', cmap(vv,:), 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', vidList(vv));
end
xlabel('time (min)'); ylabel('mean DQ-BSA object area (px)')
title([genoTypePrefix, ' droplet size over chase'])
legend('Location','eastoutside','Interpreter','none'); box off
%% object count per frame, segmentation sanity check
figure; hold on
for vv=1:length(vidList)
    thisVid=meanTable(meanTable.video==vidList(vv), :);
    plot(thisVid.timeMin, thisVid.nObj, '-o', 'Color', cmap(vv,:), 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', vidList(vv));
end
xlabel('time (min)'); ylabel('n DQ-BSA objects detected')
title('segmentation stability'); legend('Location','eastoutside','Interpreter','none'); box off