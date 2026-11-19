%Pulse Chase Pipeline 




%% load videos and set parameters




%% initializing relevant variables and loading video
addpath('/Users/pennerc/Documents/Documents - BGS-STU23-138/MATLAB/GPNMB_Global_Analysis/GPNMB-Global-Analysis');
% first upload relevant video

homePath= '/Volumes/PC60/iMicroVids/PulseChase/7_16_26/KO/';


genoTypePrefix='KO';
imageTypePrefix='pulseChase';

outputSuffix='7_6_26';

frameRate=5; %one image every 5 minutes



%Mask Making factors
minSize=125; % 125 for 10x 250 for 10x lmao this is the minimum size an object can be (some dapi debris always present)
maxSize=700;   %original run 600 % sometimes multiple cells become clumped (as in a swarm) I just remove these
circCut=0.33; % originally did .05 this is a cutoff for the minimum circularity of a detected object, less circular objects tend to be debris n whatnot.  
detectSensitivity=.5; % for now this is meaningless insofar as we are using OTSU's this is the sensitivity of our initial binarization step. I've found .6 is a good middle ground but decrease if there are a lot of artifactual objects
strelRadiusSize=3; %this is the pixel radius of our structuring element. We first apply an erosion of all pixels to remove noise and then proceed with dilation procedure
%a larger strel Radius Size will cause objects that are further apart ( 4
%vs 5 vs 6 pixels) to become aggregated
demoPlot=false; %this will just create a figure with the masks at each of our cleanup stages





  %Bodipy Staining
    
    bodipMaxSize=30000; % 
    bodipMinSize=400; %  edges of neurites sometimes separated by strel
    bodipDetectSensitivity=.5;
    bodipRadiusSize=5; %here strel type is always arbitrary
    bodipCircCut=0;
    bodipStrelType='arbitrary';


% % trajectory filtering features
% maxJump=30; %this is the upper limit of how much a nucleus can move in one frame
% confThreshold=.05; %this is the lower limit of our confidence in tracking a nucleus from one point to another.
% areaThreshold=1; %this is how much of a difference in size from frame to frame you'll tolerate this will be protection against fusion or splitting events
% 

% %Scratch analysis factors
% % Set threshold in pixels ( how far from scratch counts as "near")
% proximityThreshold = 5;
% % genoTypePrefix='KO';
% % imageTypePrefix='Scratch';


allFilesAtPlay = dir(fullfile(homePath, '*.mp4'));
allFilesAtPlay = allFilesAtPlay(~startsWith({allFilesAtPlay.name}, '._')); % mac keeps these weird ghost files with ._ prefix

allFilesAtPlay= allFilesAtPlay(  contains( {allFilesAtPlay.name}, genoTypePrefix)  & ...
(contains( {allFilesAtPlay.name}, imageTypePrefix)))




nucMaskCollect=cell(1,length(allNames)); nucFeatureCollect= cell(1,length(allNames));rawNucCollect= cell(1,length(allNames));
bodipMaskCollect=cell(1,length(allNames)); bodipFeatureCollect= cell(1,length(allNames)); bodipRawCollect= cell(1,length(allNames));



%% start of for Loop through all vids

 for vidNum=1: length(allFilesAtPlay)            
        
        vidName=[homePath,'/', allFilesAtPlay(vidNum).name];
 
         %KO_Scratch_aB (1).mp4
         vidObj = VideoReader(vidName);
    

    
        % loading objects into frames
        
        disp('loading video data, this will take a few seconds')
        
        %preallocating objects
        
        numFrames = floor(vidObj.Duration * vidObj.FrameRate);
        dapiFrames = cell(numFrames, 1);
        bodipFrames= cell(numFrames,1)
        refFrames=cell(numFrames, 1);
        frameIdx = 1;

        while hasFrame(vidObj)
            rgbFrame = readFrame(vidObj);
        
            % Extract the blue channel (DAPI)
            blueChannel = rgbFrame(:,:,3);  % Blue is the 3rd channel in RGB


            dapiFrames{frameIdx} = blueChannel;
            refFrames{frameIdx}=rgbFrame;



        
            frameIdx = frameIdx + 1;
        end
        



        
        maskCollect=cell(numFrames, 1); %we will collect both masks and info about them
        frameDataCollect=cell(numFrames, 1);
        objNumCollect=nan(1,numFrames);
        
        
        
            if ~demoPlot
                    frames2Assess=length(dapiFrames);
                    disp('video loaded, generating masks for each frame of movie')
        
            elseif demoPlot
                frames2Assess=6;
                disp('outputting image of binarized masks')
                figure
        
            end
        
        
                for bb=1: frames2Assess
                    if ~isempty(dapiFrames{bb})
                [maskCollect{bb},CC,frameDataCollect{bb}] = pulseChaseMaskMakerDapi(dapiFrames{bb},minSize, maxSize, circCut, detectSensitivity,strelRadiusSize,demoPlot);
                objNumCollect(bb)=height(frameDataCollect{bb});
                   
                    end

                end
        







%% mask making (anchor only to lipid)






%% analysis, of changes over time, normalize to baseline

% secondary quesitons include morphology changes and single cell tracking 






%% final use glme output is change in fluoro*time
