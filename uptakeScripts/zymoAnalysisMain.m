%zymo uptake main

%% 1) loading everything and setting initial parameters
%load images (these are spearated into three channels, C1 Trem2 C2 aB C3
%Nuclei)
%this is designed for a 20x field of view. 

addpath('/Users/pennerc/Documents/Documents - BGS-STU23-138/MATLAB/GPNMB_Global_Analysis/GPNMB-Global-Analysis');


%Link the folder containing all image 

  

%this is obvi not elegant but I think of it like having to pull three
%separate switches each run
set1=true;
set2=false;


if set1
homePath='/Volumes/PC60/phagocytosisAssays/Zymo/Set1/RowD_Zymo';
elseif set2
homePath='/Volumes/PC60/phagocytosisAssays/Zymo/Set2/ZYMO'; 
end

%'/Volumes/PC60/ABeta_Uptake/files_set2';
% '/Volumes/PC60/ABeta_Uptake/viewable_tiffs_07302025'



% we'll iterate through the directory
allFilesAtPlay = dir(fullfile(homePath, '*.tif'));
allFilesAtPlay = allFilesAtPlay(~startsWith({allFilesAtPlay.name}, '._')); % mac keeps these weird ghost files with ._ prefix



%these values will change a bit between sets due to quality of images, 

% mask making factors


    nucMinSize=250; %Dapi Debris
    nucMaxSize=8000;
    nuccircCut=0.5;% this is a cutoff for the minimum circularity of a detected object, less circular objects tend to be debris n whatnot. 
    nucStrelType='disk';
    nucStrelRadiusSize=3; %this is the pixel radius of our structuring element. We first apply an erosion of all pixels to remove noise and then proceed with dilation procedure



    %we can use the same pipeline for aB just shifting size threshold and
    %removing our Circularity requirement
    zymoCircCut=0; % this is a cutoff for the minimum circularity of a detected object, less circular objects tend to be debris n whatnot.  
    zymoMaxSize=50000; % 
    zymoMinSize=50;
    zymoStrelType='arbitrary'; %don't want to bias shape of aB 
    zymoStrelRadius=3;
    zymoSensitivity=.5;



    %the Trem2 stain is a bit fussier, originally I did a mixed approach
    %but actually just local separation with super high intensity works
    %well
    
    T2maxSize=30000; % 
    T2minSize=400; %  edges of neurites sometimes separated by strel
    T2detectSensitivity=.5;
    T2strelRadiusSize=5; %here strel type is always arbitrary
    T2CircCut=0;
    T2strelType='arbitrary';


%plotting parameters
    demoPlot=false; %for initial testing of factors this will output masks for each of the fluorescent channels  to facilitate tweaking
    totPlot=false; % for outputting the final overlayed mask in comparison to the originial image. Imp to note this is pre fluorescent image. 

% analysis parameters
   aBDistCutoff=70; % this is the allowable distance (in pixels) for an aB object from a given nucleus. This was generated based on visual assesment
   %of distance distribution (calculated in section 4 below)
   zymoFluoroCutoff=75; %this is the minimum fluoro for an object to be considered "real" It was calaculated through visual assesment of the fluoro distribution
   % this was output in section 3 below




% outputting the image names to iterate through

nameCollect=[];
for bb=1:length(allFilesAtPlay)
    nameAtPlay=allFilesAtPlay(bb).name;
    dashPoint= strfind(nameAtPlay, '_');  
    nameCollect{bb}=nameAtPlay(1:dashPoint(end)-1);

end

   
allNames=unique(nameCollect);

%holding variables for our mask feature outputs
nucMaskCollect=cell(1,length(allNames)); nucFeatureCollect= cell(1,length(allNames));rawNucCollect= cell(1,length(allNames));
zymoMaskCollect=cell(1,length(allNames)); zymoFeatureCollect= cell(1,length(allNames)); rawZymoCollect= cell(1,length(allNames));
T2MaskCollect=cell(1,length(allNames)); T2FeatureCollect= cell(1,length(allNames)); rawT2Collect= cell(1,length(allNames));




%% 2) capturing basic masks for all channels as well as raw data
    %iterating through each name. this will contain three different
    %channels this will take a while



%initializing wait bar
f = waitbar(0, 'Starting');
  
 for dd=82:96%length(allNames)
         %I'm fairly certain I could just iterate through in sets of 3, but
         %I want to be agnostic to ordering in case that's not always true.
     nameBaseAtPlay=allNames{dd};


     if set1
     nuclearImAtPlay=imread( [homePath,'/',nameBaseAtPlay,'_DRAQ5.tif']); rawNucCollect{dd}=nuclearImAtPlay; % I grab the raw images for downsttream plotting and intensity filtering
      Trem2ImAtPlay=imread( [homePath,'/',nameBaseAtPlay,'_TREM2.tif']);   rawT2Collect{dd}=Trem2ImAtPlay;

     zymoImAtPlay=imread( [homePath,'/',nameBaseAtPlay,'_ZYMO.tif']);      rawZymoCollect{dd}=zymoImAtPlay; %weird misname thing


    elseif set2
      nuclearImAtPlay=imread( [homePath,'/',nameBaseAtPlay,'_DRAQ5.tif']); rawNucCollect{dd}=nuclearImAtPlay; % I grab the raw images for downsttream plotting and intensity filtering
      Trem2ImAtPlay=imread( [homePath,'/',nameBaseAtPlay,'_TREM2.tif']);   rawT2Collect{dd}=Trem2ImAtPlay;

     zymoImAtPlay=imread( [homePath,'/',nameBaseAtPlay,'_LYSO.tif']);      rawZymoCollect{dd}=zymoImAtPlay; %weird misname thing


     end


    [nucMaskCollect{dd},~,nucFeatureCollect{dd}] = aBNuclearMaskMaker(nuclearImAtPlay,nucMinSize, nucMaxSize, nuccircCut,nucStrelRadiusSize,nucStrelType,demoPlot);

    
if demoPlot
  title('nuclear')  
end

    [zymoMaskCollect{dd},~,zymoFeatureCollect{dd}] = zymoMaskMaker(zymoImAtPlay,zymoMinSize, zymoMaxSize, zymoCircCut,zymoStrelRadius,zymoStrelType,zymoSensitivity, demoPlot);

if demoPlot
      title('aB')  

end
    [T2MaskCollect{dd},~,T2FeatureCollect{dd}] = aBNuclearMaskMakerTrem2Specific(Trem2ImAtPlay,T2minSize, T2maxSize, T2CircCut,T2strelRadiusSize,T2strelType, T2detectSensitivity, demoPlot);

if demoPlot

    title('trem2')  
    
end

    waitbar(dd/length(allNames), f, sprintf('Progress: %d %%', floor( (dd/length(allNames))  *100)));
    pause(0.1);



     totIm = cat(3, zymoImAtPlay, Trem2ImAtPlay, nuclearImAtPlay);     % 3D array HxWx3
     maskImTot = cat(3, zymoMaskCollect{dd}, T2MaskCollect{dd}, nucMaskCollect{dd});  
    if totPlot  
     figure
     imshow(totIm)
     % imshowpair(totIm,maskImTot,'montage')
     % title(["over view of  ", allNames{dd}])
    end



 end
 disp('nice, everything loaded!')
close(f)
        % so this seems to mostly work.... because I'm using otsu for some of these it can pick up very low signal. So next we'll create an array of mean fluorescence for our TREM2
   %and aBeta signals, this will help us create a fluorescence filter. I actually think a lot of these low signals, esp the TREm2 are real.  



%% 3) loading all Fluoresences to generate a cutoff (you can do this automatically with just some kind of std metric but I think it's better to finesse it)


aSynFluoroCollect=[];
T2FluoroCollect=[];
            for tt=1:length(allNames)
                    rawAsynIm= rawZymoCollect{tt};
                    rawT2Im=rawT2Collect{tt};
                    zymoFeaturesAtPlay=zymoFeatureCollect{tt};
                    T2FeaturesAtPlay=T2FeatureCollect{tt};

                    for zz=1:height(zymoFeaturesAtPlay)
                        aSynFluoroCollect(1+length(aSynFluoroCollect))= nanmean(rawAsynIm(zymoFeaturesAtPlay.PixelIdxList{zz})  ); %simply grab the pixel values found in the space of our mask

                    end

                    for zz=1:height(T2FeaturesAtPlay)

                        T2FluoroCollect(1+length(T2FluoroCollect))= nanmean(rawT2Im(T2FeaturesAtPlay.PixelIdxList{zz})  ); %likewise for trem2

                    end
            end

% based on just looking it seems like a Fluoro Cut off of 10 is reasonable
% for both channels. There is a clear high degree of low signal Trem2, less
% so for aB. 

             % 
             % figure
             % histogram(T2FluoroCollect)

            figure
            histogram(aSynFluoroCollect)

hold on


           




            % 
            %  figure
            % histogram(T2FluoroCollect(T2FluoroCollect<25))




%% 4) now doing a similar analysis to establish analysis parameters
% by looking at the distance from nucleus to putative aB object IMO this
% didn't really work, I think this really underestimated the distance to
% real objects. I basically generated this value by looking at a lot of
% different outputs across various conditions. 
% 
% 
% nucDistanceCollect=[];
%             for tt=1: length(allNames)
% 
%                     rawABIm= rawABCollect{tt};
%                     aBFeaturesAtPlay=aBFeatureCollect{tt};
%                     abMaskAtPlay=aBMaskCollect{tt};
%                     nucFeaturesAtPlay=nucFeatureCollect{tt};
% 
% 
%                        if ~isempty(aBFeaturesAtPlay)
%                         %eucDist=pdist2(nucFeaturesAtPlay.Centroid, aBFeaturesAtPlay.Centroid  );
%                         distMap = bwdist(abMaskAtPlay, 'euclidean');
% 
%                             rows = round(nucFeaturesAtPlay.Centroid(:,2));  % y = row
%                             cols = round(nucFeaturesAtPlay.Centroid(:,1));  % x = col
% 
%                             idx = sub2ind(size(distMap), rows, cols);
%                             eucDist = distMap(idx);
%                             eucDist = eucDist(:).'; %so this isn't perfect but I'm just using this for my rough cut offs (ie what I'll consider within a reasonable distance of a nucleus
%                             %We take the distance from the nuclei centroid
%                             %to the nearest positive location in the aB
%                             %mask.
% 
% 
%                             nucDistanceCollect=[nucDistanceCollect,eucDist];
%                        else
%                             nucDistanceCollect= [nucDistanceCollect, nan(1, height(nucFeaturesAtPlay))];
%                        end                        
%             end
% 
% 
% 
% 
%             figure
% 
%             histogram(nucDistanceCollect(nucDistanceCollect<100));
% 
% 

            
            %  seems like 1 or 2 is the right cutoff really most seem to be
            %  directly adjacent to an aB 


            %% 5) alright analysis time finally





        zymoSpotsCellProx=nan(1, length(allNames)); % this will be the number of asyn spots within radius
        avgZymoGluoroCellProx=nan(1, length(allNames)); %likewise for aB fluro
        avgZymoAreaCellProx=nan(1, length(allNames)); %likewise for aB area

        numZymoCellProx= nan(1, length(allNames)); 

        avgZymoFluoroTot=nan(1, length(allNames)); 
        avgZymoAreaTot=nan(1, length(allNames)); 
        numZymoSpotsTot=nan(1, length(allNames)); 

        numNuclei=nan(1, length(allNames));



 
%initializing wait bar
f = waitbar(0, 'Starting');


            for tt=1: length(allNames)

                nucFeaturesAtPlay=nucFeatureCollect{tt}; % loop through nuclei in a given image
    
                rawAsynIm= rawZymoCollect{tt}; % to grab fluoresence 
                zymoFeaturesAtPlay=zymoFeatureCollect{tt};% for size and location
                zymoMaskAtPlay=zymoMaskCollect{tt}; % for mapping
                aBDistMap = bwdist(zymoMaskAtPlay, 'euclidean');
    
    
                 rawTremIm= rawT2Collect{tt}; % to grab fluoresence 
                 tremFeaturesAtPlay=T2FeatureCollect{tt};% for size and location
                 tremMaskAtPlay=T2MaskCollect{tt};
                 tremDistMap = bwdist(tremMaskAtPlay, 'euclidean');


                zymoAreaCollectNuc=nan(1,height(nucFeaturesAtPlay));
                zymoFluoroCollectNuc=nan(1,height(nucFeaturesAtPlay));
                zymoTrueFalseNuc=false(1,height(nucFeaturesAtPlay));
                zymoCountCollect=false(1,height(nucFeaturesAtPlay));


                zymoAreaCollectTot=nan(1,height(nucFeaturesAtPlay));
                zymoFlouroCollectTot=nan(1,height(nucFeaturesAtPlay));
                zymoTrueFalseTot=false(1,height(nucFeaturesAtPlay));


                T2AreaCollectNucAB=nan(1,height(nucFeaturesAtPlay));
                T2fluoroCollectNucAB=nan(1,height(nucFeaturesAtPlay));
                T2DistCollectNucAB=nan(1,height(nucFeaturesAtPlay));

                T2AreaCollectNucNonAB=nan(1,height(nucFeaturesAtPlay));
                T2fluoroCollectNucNonAB=nan(1,height(nucFeaturesAtPlay));
                T2DistCollectNucNonAB=nan(1,height(nucFeaturesAtPlay));


                blankMap=false(size(rawAsynIm));

%here we are outputting the distance of each pixel on a map from the
%nearest aB object or T2 object
                                        mapHoldZymo=cell(1, height(zymoFeaturesAtPlay));
                                     for cc=1:height(zymoFeaturesAtPlay)
                                       mapperSupreme= blankMap; mapperSupreme(zymoFeaturesAtPlay.PixelIdxList{cc})=true;
                                       mapHoldZymo{cc}=bwdist(mapperSupreme, 'euclidean');
                                     end


                                             mapHoldT2=cell(1, height(tremFeaturesAtPlay));
                                     for cc=1:height(tremFeaturesAtPlay)
                                       mapperSupreme= blankMap; mapperSupreme(tremFeaturesAtPlay.PixelIdxList{cc})=true;
                                       mapHoldT2{cc}=bwdist(mapperSupreme, 'euclidean');
                                     end


               

                                %first grab all aB objects within 2 pixels of a given nucleus (anchored to cells) 
                      for dd=1:height(nucFeaturesAtPlay)
  
% nucleus should always be in the midst of a trem2 object, unless it
% doesn't exist

                                      rows = round(nucFeaturesAtPlay.Centroid(dd,2));  % y = row
                                      cols = round(nucFeaturesAtPlay.Centroid(dd,1));  % x = col
                                      idx = sub2ind(size(aBDistMap), rows, cols);

% this is a little dangerous since this will always output SOMETHING but we
% can filter based on distance and fluoresence downstream
                               
                                     eucDistTrem=[];
                                     for gg=1:height(tremFeaturesAtPlay)
                                       eucDistTrem(gg)=mapHoldT2{gg}(idx);
                                     end
if ~isempty(tremFeaturesAtPlay)
%trem2DistRough=pdist2(nucFeaturesAtPlay.Centroid(dd,:), tremFeaturesAtPlay.Centroid  ); % get the distance from our queried nucleus to each trem2object
                                [~,trem2Idx]= min(eucDistTrem);
                                T2FluoroAtPlay=nanmean(rawTremIm(tremFeaturesAtPlay.PixelIdxList{trem2Idx}));
                                T2AreaAtPlay=tremFeaturesAtPlay.Area(trem2Idx);
                                T2DistAtPlay=eucDistTrem(trem2Idx);

                                if T2DistAtPlay >25 || T2FluoroAtPlay <10 % really hard for me to pick a fluoresence threshold, this is weirdly homogenous? Maybe some kind of issue
                                    T2DistAtPlay=nan;
                                    T2AreaAtPlay=nan;
                                    T2FluoroAtPlay=nan;
                                end


elseif isempty(tremFeaturesAtPlay)

                                    T2DistAtPlay=nan;
                                    T2AreaAtPlay=nan;
                                    T2FluoroAtPlay=nan;

end

                           if ~isempty(zymoFeaturesAtPlay) % if there are any aB objects to begin with





                                        eucDistZymo=[];
                                     for gg=1:height(zymoFeaturesAtPlay)
                                       eucDistZymo(gg)=mapHoldZymo{gg}(idx);
                                     end

                                    



                                    %step 5A making maps for each aB object... probably an easier way  to solve this

                                     abIdx=find(eucDistZymo'<=aBDistCutoff)  ; % first we apply the distance cutoff

                
                               % eucDist=pdist2(nucFeaturesAtPlay.Centroid(dd,:), aBFeaturesAtPlay.Centroid  ); % get the distance from our queried nucleus to each aBObject
                
                                         if ~isempty(abIdx) % if we have any objects that are within our distance
        
                                             aBPowerCollectLocal=nan(1, length(abIdx));
                                             aBAreaCollectLocal=nan(1, length(abIdx));
                                                %step 5B iterating through
                                                %the idx of possible aB
                                                %candidates and outputing
                                                %the area and fluoro of
                                                %each
    
                                                 for bb=1:length(abIdx) 
                                                     aBPowerCollectLocal(bb)= nanmean(rawAsynIm(zymoFeaturesAtPlay.PixelIdxList{abIdx(bb)}));
                                                     aBAreaCollectLocal(bb)= ((zymoFeaturesAtPlay.Area(abIdx(bb))));
                                                 end
                                                 % aBAreaCollectLocal(aBPowerCollectLocal <=10)=nan;
                                                 % [~,candidateAbObjectIdx]= max(aBAreaCollectLocal);
            
                                                     if sum(aBPowerCollectLocal >=zymoFluoroCutoff)>0 % if there are objects above our threshold
                                                        zymoAreaCollectNuc(dd)=nanmean(aBAreaCollectLocal(aBPowerCollectLocal>=zymoFluoroCutoff));
                                                        zymoFluoroCollectNuc(dd)=nanmean(aBPowerCollectLocal(aBPowerCollectLocal>=zymoFluoroCutoff));
                                                        zymoCountCollect(dd)=length(abIdx);

                                                        zymoTrueFalseNuc(dd)=true;


                                                        T2fluoroCollectNucAB(dd)= T2FluoroAtPlay;
                                                        T2AreaCollectNucAB(dd)=T2AreaAtPlay;
                                                        T2DistCollectNucAB(dd)=T2DistAtPlay;

                                                        T2fluoroCollectNucNonAB(dd)= nan;
                                                        T2AreaCollectNucNonAB(dd)=nan;
                                                        T2DistCollectNucNonAB(dd)=nan;
                               

                                                     else
                                                        zymoAreaCollectNuc(dd)=nan; % if there are no nearby aB objects it gets zerod
                                                        zymoFluoroCollectNuc(dd)=nan;
                                                        zymoCountCollect(dd)=0;
                                                        zymoTrueFalseNuc(dd)=false;

                                                        T2fluoroCollectNucAB(dd)= nan;
                                                        T2AreaCollectNucAB(dd)=nan;
                                                        T2DistCollectNucAB(dd)=nan;

                                                        T2fluoroCollectNucNonAB(dd)= T2FluoroAtPlay;
                                                        T2AreaCollectNucNonAB(dd)=T2AreaAtPlay;
                                                        T2DistCollectNucNonAB(dd)=T2DistAtPlay;
                                                     end
                                         else
                                            zymoAreaCollectNuc(dd)=nan;
                                            zymoFluoroCollectNuc(dd)=nan; % if there are no abeta objects that pass the dist threshold everything is zeroed!
                                            zymoCountCollect(dd)=0;
                                            zymoTrueFalseNuc(dd)=false;

                                            T2fluoroCollectNucAB(dd)= nan;
                                            T2AreaCollectNucAB(dd)=nan;
                                            T2DistCollectNucAB(dd)=nan;

                                            T2fluoroCollectNucNonAB(dd)= T2FluoroAtPlay;
                                            T2AreaCollectNucNonAB(dd)=T2AreaAtPlay;
                                            T2DistCollectNucNonAB(dd)=T2DistAtPlay;
                                         end
                           else
                                 zymoAreaCollectNuc(dd)=nan;
                                zymoFluoroCollectNuc(dd)=nan; % if there are no abeta objects period everything is zeroed!
                                zymoCountCollect(dd)=0;

                                zymoTrueFalseNuc(dd)=false;

                                T2fluoroCollectNucAB(dd)= nan;
                                T2AreaCollectNucAB(dd)=nan;
                                T2DistCollectNucAB(dd)=nan;
    
                                T2fluoroCollectNucNonAB(dd)= T2FluoroAtPlay;
                                T2AreaCollectNucNonAB(dd)=T2AreaAtPlay;
                                T2DistCollectNucNonAB(dd)=T2DistAtPlay;

                           end 



% if testTime
% 
% figure
%  mapperSupreme= blankMap; mapperSupreme(nucFeaturesAtPlay.PixelIdxList{dd})=true;
% maskImTot = cat(3, mapperSupreme, tremMaskAtPlay, nan(size(tremMaskAtPlay)));  
% 
% imshowpair(maskImTot,rawTremIm,'montage')
% 
% title(['nucleus ', num2str(dd), ' T2Dist=', num2str(T2DistAtPlay), ' T2Fluoro=',num2str(T2FluoroAtPlay) ])
% end
% 
% 
% b=1;



                      end




                        for cc=1:height(zymoFeaturesAtPlay)
                      zymoFlouroCollectTot(cc)=nanmean(rawAsynIm(zymoFeaturesAtPlay.PixelIdxList{cc}));
                      zymoAreaCollectTot(cc)=zymoFeaturesAtPlay.Area(cc);
                        end


        zymoSpotsCellProx(tt)=nanmean(zymoCountCollect); % this will be the number of aB spots within our radius of a nuclei
        avgZymoGluoroCellProx(tt)=nanmean(zymoFluoroCollectNuc); %likewise for aB fluro
        avgZymoAreaCellProx(tt)=nanmean(zymoAreaCollectNuc); %likewise for aB area

        numZymoCellProx(tt)= sum(zymoTrueFalseNuc)/  (length(zymoTrueFalseNuc)); % this will be the ratio of cells with at least one aB object within their radius 

        avgZymoFluoroTot(tt)=nanmean(zymoFlouroCollectTot(zymoFlouroCollectTot>= zymoFluoroCutoff )); % this will be the number of aB spots in an image (normalized to nuclei count)
        avgZymoAreaTot(tt)=nanmean(zymoAreaCollectTot(zymoFlouroCollectTot>= zymoFluoroCutoff ) ); % this will be the average fluoresence across the whole image
        numZymoSpotsTot(tt)=  (sum(zymoFlouroCollectTot>= zymoFluoroCutoff )   / height(nucFeaturesAtPlay)  ); % likewise with area
        numNuclei(tt)=height(nucFeaturesAtPlay);

                T2AreaCollectTotAB(tt)=nanmean(T2AreaCollectNucAB);
                T2fluoroCollectAB(tt)=nanmean(T2fluoroCollectNucAB);
                T2DistCollectAB(tt)=nanmean(T2DistCollectNucAB);



          

                T2AreaCollectTotNonAB(tt)=nanmean(T2AreaCollectNucNonAB);
                T2fluoroCollectNonAB(tt)=nanmean(T2fluoroCollectNucNonAB);
                T2DistCollectNonAB(tt)=nanmean(T2DistCollectNucNonAB);
               


   waitbar(tt/length(allNames), f, sprintf('Analysis Progression: %d %%', floor( (tt/length(allNames))  *100)));
    pause(0.1);

            end



%%

    if set1 | set2
    
   
    WT= contains(allNames,'WT');
    HET=contains(allNames,'HET');
    ECD=contains(allNames,'ECD');
    KO=contains(allNames,'KO') & ~ECD;


    end







genoArray=cell(1,length(allNames));

genoArray(WT)={'WT'}; genoArray(HET)={'HET'}; genoArray(KO)={'KO'}; genoArray(ECD)={'ECD'};




% Build table for GLME
totTable = table(...
    zymoSpotsCellProx(:), ...          % Response: avg fluorescence
    avgZymoGluoroCellProx(:), ...              % Fixed effect: genotype
    avgZymoAreaCellProx(:),...
    numZymoCellProx(:), ...            % num prox cels
    avgZymoFluoroTot(:), ...             
    avgZymoAreaTot(:),...
    numZymoSpotsTot(:), ...            
    categorical(allNames(:)), ...           % Random effect grouping
    categorical(genoArray(:)),...
    numNuclei(:),...
    'VariableNames', {'zymoSpotsProx', 'avgZymoFluoroCellProx', 'avgAreaZymoSpotsCellProx','numZymoCellProx',...
    'avgZymoFluoroTot','avgZymoAreaTot','numZymoSpotsTot','randArray','Genotype','numNuclei'} ...
);

% Set WT as reference category
totTable.Genotype = reordercats(totTable.Genotype, {'WT', 'HET', 'KO', 'ECD'});


writetable(totTable, [homePath,'/zymoUptakeQuants.csv'])



%%


val2Test=zymoSpotsCellProx;
val2Test(isinf(val2Test))=nan;

val2Test(val2Test >    (nanmean(val2Test)  +4*nanstd(val2Test)  )     )=nan;

val2Test(numNuclei<20)=nan;

ylablr=' Fold Change in number of cellular aB objects';

titlr=' Fold Change in number of cellular aB objects (treatment/PBS)  ';



figure
hold on

b=bar(1,nanmean(val2Test(WT )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [.8 .2 .5]; 
ylabel(ylablr,'FontSize',20)
aSynMaskAtPlay=gca; aSynMaskAtPlay.XTickLabel=[];



b=bar(3,nanmean(val2Test(HET )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0 0.7 .25];


b=bar(5,nanmean(val2Test(KO )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0.3 0.1 .6];


b=bar(7,nanmean(val2Test((ECD) )   ))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [.9 .1 .1]; 



hold on
aSynMaskAtPlay=scatter(rand(1, sum((WT )))+.5, val2Test(WT ), 'Marker', 'o' );
aSynMaskAtPlay.CData=[.7 .3 .6]; 

b=scatter(rand(1, sum(HET ))+2.5, val2Test(HET ), 'Marker', 'o' );
b.CData(1,:) = [0 0.7 .25];

c=scatter(rand(1, sum(KO ))+4.5, val2Test(KO), 'Marker', 'o' );
c.CData(1,:) = [0.3 0.1 .6];

aSynMaskAtPlay=scatter(rand(1, sum((ECD )))+6.5, val2Test(ECD), 'Marker', 'o' );
aSynMaskAtPlay.CData=[.9 .1 .1]; 



legend({'WT','HET','KO','KO+ECD'})



title(titlr,'FontSize',30)

% statistical testing




fitglme(totTable, 'numZymoCellProx~ 1+ Genotype + (1|randArray)')


%anovan(val2Test(uM1), {genoArray(uM1)}  )





% [p, tbl, stats] = kruskalwallis(val2Test(uM1), genoArray(uM1));


%%

                
                figure
                histogram(T2FluoroCollectTot)




                T2AreaCollectTotAB(tt)=nanmean(T2AreaCollectNucAB);
                T2fluoroCollectAB(tt)=nanmean(T2fluoroCollectNucAB);
                T2DistCollectAB(tt)=nanmean(T2DistCollectNucAB);

                T2AreaCollectTotNonAB(tt)=nanmean(T2AreaCollectNucNonAB);
                T2fluoroCollectNonAB=nanmean(T2fluoroCollectNucNonAB);
                T2DistCollectNonAB=nanmean(T2DistCollectNucNonAB);



%%


%% outputting everything as a table

% Define column names
colNames = { ...
    "ImageCode", ...
    "Ratio_of_Cells_Pos_for_aB", "Ratio_of_Cells_Pos_for_aB_normalized", ...
    "avg_aB_Fluoro_within_Cells","avg_aB_Fluoro_within_Cells_normalized", ...
    "avg_aB_area_within_Cells","avg_aB_area_within_Cells_normalized", ...
    "total_aB_objects_normalized_to_nuclei_count","total_aB_objects_normalized_to_nuclei_count relative to PBS", ...
    "avg_aB_Fluoro_Total_Image","avg_aB_Fluoro_Total_Image_normalized", ...
    "avg_aB_Area_Total_Image","avg_aB_Area_Total_Image_normalized", ...
    "Trem2_area_Pos_aB_Cells", ...
    "Trem2_area_Neg_aB_Cells", ...
    "Trem2_Fluoro_Pos_aB_Cells", ...
    "Trem2_Fluoro_Neg_aB_Cells", ...
    "Treatment",...
    "Genotype"
};


colNames = [colNames{:}];


% Build table
allData = table( ...
    allNames(:), ...
    numZymoCellProx(:), numAbCellProxNorm(:), ...
    avgZymoGluoroCellProx(:),avgAbFluoroCellProxNorm(:), ...
    avgZymoAreaCellProx(:),avgAreaAbSpotsCellProxNorm(:), ...
    numZymoSpotsTot(:),numAbSpotsTotNorm(:), ...
    avgZymoFluoroTot(:),avgAbFluoroTotNorm(:), ...
    avgZymoAreaTot(:), avgAreaAbTotNorm(:),...
    T2AreaCollectTotAB(:), ...
    T2AreaCollectTotNonAB(:), ...
    T2fluoroCollectAB(:), ...
    T2fluoroCollectNonAB(:), ...
    treatArray(:),...
    genoArray(:),...
    'VariableNames', colNames ...
);


% Save to file


    if set1
    outputPath = '/Users/pennerc/Desktop/allDataResultsSet1.xlsx';
    writetable(allData, outputPath);
    elseif set2
    outputPath = '/Users/pennerc/Desktop/allDataResultsSet2.xlsx';
    writetable(allData, outputPath);
    elseif set3
    outputPath = '/Users/pennerc/Desktop/allDataResultsSet3.xlsx';
    writetable(allData, outputPath);
    end


%%

set1Data = readtable('/Volumes/PC60/phagocytosisAssays/Zymo/Set1/RowD_Zymo/zymoUptakeQuants._set1csv.csv');
set2Data = readtable('/Volumes/PC60/phagocytosisAssays/Zymo/Set2/ZYMO/zymoUptakeQuants_set2.csv');



val2TestSet1 = set1Data.numZymoSpotsTot; set1Geno = set1Data.Genotype;


normVal2TestSet1 = val2TestSet1 / nanmean(val2TestSet1(strcmp(set1Geno,'WT')));

% normVal2TestSet1( normVal2TestSet1 > nanmean(normVal2TestSet1)+ 4* nanstd(normVal2TestSet1) )=nan;


val2TestSet2 = set2Data.numZymoSpotsTot; set2Geno = set2Data.Genotype;
normVal2TestSet2 = val2TestSet2 / nanmean(val2TestSet2(strcmp(set2Geno,'WT')));

% normVal2TestSet2( normVal2TestSet2 > nanmean(normVal2TestSet2)+ 4* nanstd(normVal2TestSet2) )=nan;


val2TestSet1(val2TestSet1> (nanmean(val2TestSet1)+4*std(val2TestSet1)))=nan;
val2TestSet2(val2TestSet2> (nanmean(val2TestSet2)+4*std(val2TestSet2)))=nan;



normVal2TestSet1(normVal2TestSet1> (nanmean(normVal2TestSet1)+4*std(normVal2TestSet1)))=nan;
normVal2TestSet2(normVal2TestSet2> (nanmean(normVal2TestSet2)+4*std(normVal2TestSet2)))=nan;


ylablr='Normalized Zymosan Fluorescence Across Genotypes';
titlr='Normalized Zymosan fluorescence';



set1Data.normVal = val2TestSet1%normVal2TestSet1;
set2Data.normVal = val2TestSet2%normVal2TestSet2;

% Add Set label
set1Data.Set = repmat("Set1", height(set1Data), 1);
set2Data.Set = repmat("Set2", height(set2Data), 1);

totTable = [set1Data; set2Data];
totTable.Set = categorical(totTable.Set);

geno = totTable.Genotype; val2Test = totTable.normVal;
isSet1 = totTable.Set == "Set1";
isSet2 = totTable.Set == "Set2";

WT  = contains(geno,'WT');
ECD = contains(geno,'ECD');
HET = contains(geno,'HET');
KO  = ~contains(geno,'ECD') & ~HET & contains(geno,'KO');

% Colors per genotype
wtCol  = [.8 .2 .5];
hetCol = [0 0.7 .25];
koCol  = [0.3 0.1 .6];
ecdCol = [0 0.2 .9];

figure; hold on

% Bars
b = bar(1, nanmean(val2Test(WT)));  b.FaceColor='flat'; b.FaceAlpha=.3; b.BarWidth=1.5; b.CData(1,:)=wtCol;
b = bar(3, nanmean(val2Test(HET))); b.FaceColor='flat'; b.FaceAlpha=.3; b.BarWidth=1.5; b.CData(1,:)=hetCol;
b = bar(5, nanmean(val2Test(KO)));  b.FaceColor='flat'; b.FaceAlpha=.3; b.BarWidth=1.5; b.CData(1,:)=koCol;
b = bar(7, nanmean(val2Test(ECD))); b.FaceColor='flat'; b.FaceAlpha=.3; b.BarWidth=1.5; b.CData(1,:)=ecdCol;

ylabel(ylablr, 'FontSize', 30);
a = gca; a.XTickLabel = [];

title(titlr,'FontSize',30)

% plottin 
groups    = {WT & ~ECD, HET & ~ECD, KO & ~ECD, ECD};
xPos      = [.5, 2.5, 4.5, 6.5];
colors    = {wtCol, hetCol, koCol, ecdCol};

for gg = 1:4
    mask = groups{gg};
    col  = colors{gg};

    % Set 1 
    idx1 = mask & isSet1;
    scatter(rand(1, sum(idx1)) + xPos(gg), val2Test(idx1), ...
        "filled", 'Marker', 'o', 'CData', repmat(col, sum(idx1), 1));

    % Set 2 
    idx2 = mask & isSet2;
    scatter(rand(1, sum(idx2)) + xPos(gg), val2Test(idx2), ...
        "filled", 'Marker', 'square', 'CData', repmat(col, sum(idx2), 1));
end

% GLME 
totTable.Genotype = categorical(totTable.Genotype);
totTable.Genotype = reordercats(totTable.Genotype, ["WT","HET","KO","ECD"]);

totTable.randArray=categorical(totTable.randArray);




glme = fitglme(totTable, ...
    'numZymoSpotsTot ~ 1 + Genotype +numNuclei+ Set + (1|randArray)', ...
    'Distribution', 'Normal', ...
    'Link', 'identity');
disp(glme)


