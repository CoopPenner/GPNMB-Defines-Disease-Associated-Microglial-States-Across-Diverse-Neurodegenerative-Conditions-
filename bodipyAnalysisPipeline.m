% alright at it again, this time for lipid signal




%% 1) loading everything and setting initial parameters
%load images (these are spearated into 2 channels, BOdipy and Dapi

%this is designed for a 20x field of view. 

addpath('/Users/pennerc/Documents/Documents - BGS-STU23-138/MATLAB/GPNMB_Global_Analysis/GPNMB-Global-Analysis');


%Link the folder containing all image 

homePath= '/Volumes/PC60/lipidStuff/untouched tifs';

%'/Volumes/PC60/lipidStuff/untouched tifs';  

%/Volumes/PC60/lipidStuff/processed_8bit


% we'll iterate through the directory
allFilesAtPlay = dir(fullfile(homePath, '*.tif'));
allFilesAtPlay = allFilesAtPlay(~startsWith({allFilesAtPlay.name}, '._')); % mac keeps these weird ghost files with ._ prefix



%these values will change a bit between sets due to quality of images, 

% mask making factors


    nucMinSize=250; %Dapi Debris
    nucMaxSize=10000;
    nuccircCut=0.5;% this is a cutoff for the minimum circularity of a detected object, less circular objects tend to be debris n whatnot. 
    nucStrelType='disk';
    nucStrelRadiusSize=3; %this is the pixel radius of our structuring element. We first apply an erosion of all pixels to remove noise and then proceed with dilation procedure



    %we can use the same pipeline for aB just shifting size threshold and
    %removing our Circularity requirement
    bodipCircCut=0; % this is a cutoff for the minimum circularity of a detected object, less circular objects tend to be debris n whatnot.  
    bodipMaxSize=30000; % 
    bodipMinSize=60;
    bodipStrelType='arbitrary'; %don't want to bias shape of aB 
    bodipStrelRadiusSize=10;

bodipDistCutoff=2;
bodipFluoroCutoff=8000;



%plotting parameters
    demoPlot=true; %for initial testing of factors this will output masks for each of the fluorescent channels  to facilitate tweaking
    totPlot=true; % for outputting the final overlayed mask in comparison to the originial image. Imp to note this is pre fluorescent image. 



% outputting the image names to iterate through

nameCollect=[];
for bb=1:length(allFilesAtPlay)
    nameAtPlay=allFilesAtPlay(bb).name;
    dashLoc=find(nameAtPlay=='_');

    nameCollect{bb}=nameAtPlay(1:dashLoc-1);

end

   
allNames=unique(nameCollect);

%holding variables for our mask feature outputs
nucMaskCollect=cell(1,length(allNames)); nucFeatureCollect= cell(1,length(allNames));rawNucCollect= cell(1,length(allNames));
bodipMaskCollect=cell(1,length(allNames)); bodipFeatureCollect= cell(1,length(allNames)); rawBodipCollect= cell(1,length(allNames));


%% 2) capturing basic masks for all channels as well as raw data
    %iterating through each name. this will contain 2 different



%initializing wait bar
f = waitbar(0, 'Starting');
  
 for dd= 3%[1,2,3,5,8,9] %1:length(allNames)
         %I'm fairly certain I could just iterate through in sets of 3, but
         %I want to be agnostic to ordering in case that's not always true.
     nameBaseAtPlay=allNames{dd};

          %for set 1/2

    nuclearImAtPlay=imread( [homePath,'/',nameBaseAtPlay,'_DAPI.tif']); rawNucCollect{dd}=nuclearImAtPlay; % I grab the raw images for downsttream plotting and intensity filtering
      bodipImAtPlay=imread( [homePath,'/',nameBaseAtPlay,'_Bodipy.tif']);      rawBodipCollect{dd}=bodipImAtPlay;



    [nucMaskCollect{dd},~,nucFeatureCollect{dd}] = bodipNuclearMaskMaker(nuclearImAtPlay,nucMinSize, nucMaxSize, nuccircCut,nucStrelRadiusSize,nucStrelType, 'global', demoPlot);

    
if demoPlot
  title('nuclear')  
end

% we will grab the diffuse bodipy signal then try to evaluate puncta after.
% 
    [bodipMaskCollect{dd},~,bodipFeatureCollect{dd}] = bodipNuclearMaskMaker(bodipImAtPlay,bodipMinSize, bodipMaxSize, bodipCircCut,bodipStrelRadiusSize,bodipStrelType,'adaptive', demoPlot);

if demoPlot
      title('Bodipy')  

end

    waitbar(dd/length(allNames), f, sprintf('Progress: %d %%', floor( (dd/length(allNames))  *100)));
    pause(0.1);



     totIm = cat(3, bodipImAtPlay, nan(size(bodipImAtPlay)), nuclearImAtPlay);     % 3D array HxWx3
     maskImTot = cat(3,  nan(size(bodipMaskCollect{dd})),bodipMaskCollect{dd}, nucMaskCollect{dd});  
    if totPlot  
     figure
     imshowpair(totIm,maskImTot,'montage')
     title(["over view of  ", allNames{dd}])
    end



 end
 disp('nice, everything loaded!')
close(f)
        % so this seems to mostly work.... because I'm using otsu for some of these it can pick up very low signal. So next we'll create an array of mean fluorescence for our TREM2
   %and aBeta signals, this will help us create a fluorescence filter. I actually think a lot of these low signals, esp the TREm2 are real.  


%%
isWT=   (strcmp(allNames, 'Series2')    )| strcmp(allNames, 'Series1') |...
    strcmp(allNames, 'Series3') | strcmp(allNames, 'Series4')|...
    strcmp(allNames, 'Series25') | strcmp(allNames, 'Series26')|...
    strcmp(allNames, 'Series27') | strcmp(allNames, 'Series28')|...
    strcmp(allNames, 'Series37') | strcmp(allNames, 'Series38')|...
    strcmp(allNames, 'Series39') | strcmp(allNames, 'Series40');



bodipyFluoroCollect=[];
            for tt=  1: length(allNames)
                    rawBodipIm=rawBodipCollect{tt};
                    bodipFeaturesAtPlay=bodipFeatureCollect{tt};
    
                    for zz=1:height(bodipFeaturesAtPlay)
                        bodipyFluoroCollect(1+length(bodipyFluoroCollect))= nanmean(rawBodipIm(bodipFeaturesAtPlay.PixelIdxList{zz})  ); %simply grab the pixel values found in the space of our mask
    
                    end
    
  
            end

       

            % 
              figure
              histogram(bodipyFluoroCollect, 'BinWidth',500)
hold on
             histogram(saveyJones, 'BinWidth',500)

%%

%% 5) alright analysis time finally



            % we will iterate through each name, and then anchor analyses
            % on a per nucleus basis. 

        numBodipSpots=nan(1,length(allNames)); % this will be the number of aB spots within our radius of a nuclei
        avgBodipFluoroTot=nan(1,length(allNames)); %likewise for aB fluro
        avgBodipAreaTot=nan(1,length(allNames)); %likewise for aB area

        numBodipCellProx=nan(1,length(allNames)); % this will be the ratio of cells with at least one aB object within their radius 

        avgBodipFluroCellProx=nan(1,length(allNames)); % this will be the average fluoresence across the whole image
        avgBodipAreaCellProx=nan(1,length(allNames)); % likewise with area

       
%initializing wait bar
f = waitbar(0, 'Starting');


            for tt=1: length(allNames)

                nucFeaturesAtPlay=nucFeatureCollect{tt}; % loop through nuclei in a given image
    
                rawBodipIm= rawBodipCollect{tt}; % to grab fluoresence 
                bodipFeaturesAtPlay=bodipFeatureCollect{tt};% for size and location
                bodipMaskAtPlay=bodipMaskCollect{tt}; % for mapping
                bodipDistMap = bwdist(bodipMaskAtPlay, 'euclidean');
    
   
                bodipAreaCollectNuc=nan(1,height(nucFeaturesAtPlay));
                bodipFluoroCollectNuc=nan(1,height(nucFeaturesAtPlay));
                bodipTrueFalseNuc=false(1,height(nucFeaturesAtPlay));
                bodipCountCollect=false(1,height(nucFeaturesAtPlay));


                bodipAreaCollectTot=nan(1,height(nucFeaturesAtPlay));
                bodipFluoroCollectTot=nan(1,height(nucFeaturesAtPlay));
                bodipTrueFalseTot=false(1,height(nucFeaturesAtPlay));

                blankMap=false(size(rawBodipIm));


                                        mapHoldBodip=cell(1, height(bodipFeaturesAtPlay));
                                     for cc=1:height(bodipFeaturesAtPlay)
                                       mapperSupreme= blankMap; mapperSupreme(bodipFeaturesAtPlay.PixelIdxList{cc})=true;
                                       mapHoldBodip{cc}=bwdist(mapperSupreme, 'euclidean');
                                     end


                                %first grab all aB objects within 2 pixels of a given nucleus (anchored to cells) 
                      for dd=1:height(nucFeaturesAtPlay)
  
% nucleus should always be in the midst of a Bodipy object, unless it
% doesn't exist

                                      rows = round(nucFeaturesAtPlay.Centroid(dd,2));  % y = row
                                      cols = round(nucFeaturesAtPlay.Centroid(dd,1));  % x = col
                                      idx = sub2ind(size(bodipDistMap), rows, cols);

                           if ~isempty(bodipFeaturesAtPlay) % if there are any aB objects to begin with


                                        eucDistBodip=[];
                                     for gg=1:height(bodipFeaturesAtPlay)
                                       eucDistBodip(gg)=mapHoldBodip{gg}(idx);
                                     end

                                    %step 5A making maps for each aB object... probably an easier way  to solve this

                                     bodipIdx=find(eucDistBodip'<=bodipDistCutoff)  ; % first we apply the distance cutoff

                                
                                         if ~isempty(bodipIdx) % if we have any objects that are within our distance
        
                                             bodipPowerCollectLocal=nan(1, length(bodipIdx));
                                             bodipAreaCollectLocal=nan(1, length(bodipIdx));
                                                %step 5B iterating through
                                                %the idx of possible Bodipy
                                                %candidates and outputing
                                                %the area and fluoro of
                                                %each
    
                                                 for bb=1:length(bodipIdx) 
                                                     bodipPowerCollectLocal(bb)= nanmean(rawBodipIm(bodipFeaturesAtPlay.PixelIdxList{bodipIdx(bb)}));
                                                     bodipAreaCollectLocal(bb)= ((bodipFeaturesAtPlay.Area(bodipIdx(bb))));
                                                 end
                                        
            
                                                     if sum(bodipPowerCollectLocal >=bodipFluoroCutoff)>0 % if there are objects above our threshold
                                                        bodipAreaCollectNuc(dd)=nanmean(bodipAreaCollectLocal(bodipPowerCollectLocal>=aBFluoroCutoff));
                                                        bodipFluoroCollectNuc(dd)=nanmean(bodipPowerCollectLocal(bodipPowerCollectLocal>=aBFluoroCutoff));
                                                        bodipCountCollect(dd)=length(bodipIdx);

                                                        bodipTrueFalseNuc(dd)=true;

                                                     else
                                                        bodipAreaCollectNuc(dd)=nan; % if there are no nearby aB objects it gets zerod
                                                        bodipFluoroCollectNuc(dd)=nan;
                                                        bodipCountCollect(dd)=0;
                                                        bodipTrueFalseNuc(dd)=false;

                                                     end
                                         else
                                            bodipAreaCollectNuc(dd)=nan;
                                            bodipFluoroCollectNuc(dd)=nan; % if there are no abeta objects that pass the dist threshold everything is zeroed!
                                            bodipCountCollect(dd)=0;
                                            bodipTrueFalseNuc(dd)=false;

                                         end
                           else
                                 bodipAreaCollectNuc(dd)=nan;
                                bodipFluoroCollectNuc(dd)=nan; % if there are no abeta objects period everything is zeroed!
                                bodipCountCollect(dd)=0;

                                bodipTrueFalseNuc(dd)=false;

                           end 

                      end

                        for cc=1:height(bodipFeaturesAtPlay)
                      bodipFluoroCollectTot(cc)=nanmean(rawBodipIm(bodipFeaturesAtPlay.PixelIdxList{cc}));
                      bodipAreaCollectTot(cc)=bodipFeaturesAtPlay.Area(cc);
                        end


        avgBodipFluroCellProx(tt)=nanmean(bodipFluoroCollectNuc); %likewise for aB fluro
        avgBodipAreaCellProx(tt)=nanmean(bodipAreaCollectNuc); %likewise for aB area

        numBodipCellProx(tt)= sum(bodipTrueFalseNuc)/  (length(bodipTrueFalseNuc)); % this will be the ratio of cells with at least one aB object within their radius 

        avgBodipFluoroTot(tt)=nanmean(bodipFluoroCollectTot(bodipFluoroCollectTot>= aBFluoroCutoff )); % this will be the number of aB spots in an image (normalized to nuclei count)
        avgBodipAreaTot(tt)=nanmean(bodipAreaCollectTot(bodipFluoroCollectTot>= aBFluoroCutoff ) ); % this will be the average fluoresence across the whole image
        numBodipSpots(tt)=  (sum(bodipFluoroCollectTot>= aBFluoroCutoff )   / height(nucFeaturesAtPlay)  ); % likewise with area

         
        
   waitbar(tt/length(allNames), f, sprintf('Analysis Progression: %d %%', floor( (tt/length(allNames))  *100)));
    pause(0.1);

            end

%% setting names n whatnot

isWT=   (strcmp(allNames, 'Series2')    )| strcmp(allNames, 'Series1') |...
    strcmp(allNames, 'Series3') | strcmp(allNames, 'Series4')|...
    strcmp(allNames, 'Series13') | strcmp(allNames, 'Series14')|...
    strcmp(allNames, 'Series15') | strcmp(allNames, 'Series16')|...
    strcmp(allNames, 'Series25') | strcmp(allNames, 'Series26')|...
    strcmp(allNames, 'Series27') | strcmp(allNames, 'Series28')|...
    strcmp(allNames, 'Series37') | strcmp(allNames, 'Series38')|...
    strcmp(allNames, 'Series39') | strcmp(allNames, 'Series40');



isHet=   (strcmp(allNames, 'Series7')    )| strcmp(allNames, 'Series8') |...
    strcmp(allNames, 'Series5') | strcmp(allNames, 'Series6')|...
    strcmp(allNames, 'Series17') | strcmp(allNames, 'Series18')|...
    strcmp(allNames, 'Series19') | strcmp(allNames, 'Series20')|...
    strcmp(allNames, 'Series29') | strcmp(allNames, 'Series30')|...
    strcmp(allNames, 'Series27') | strcmp(allNames, 'Series28')|...
    strcmp(allNames, 'Series41') | strcmp(allNames, 'Series42')|...
    strcmp(allNames, 'Series43') | strcmp(allNames, 'Series44');




isKO=   (strcmp(allNames, 'Series9')    )| strcmp(allNames, 'Series10') |...
    strcmp(allNames, 'Series11') | strcmp(allNames, 'Series12')|...
    strcmp(allNames, 'Series21') | strcmp(allNames, 'Series22')|...
    strcmp(allNames, 'Series23') | strcmp(allNames, 'Series24')|...
    strcmp(allNames, 'Series33') | strcmp(allNames, 'Series34')|...
    strcmp(allNames, 'Series35') | strcmp(allNames, 'Series36')|...
    strcmp(allNames, 'Series45') | strcmp(allNames, 'Series46')|...
    strcmp(allNames, 'Series47') | strcmp(allNames, 'Series48');


isPBS=   (strcmp(allNames, 'Series1')    )| strcmp(allNames, 'Series2') |...
    strcmp(allNames, 'Series3') | strcmp(allNames, 'Series10')|...
    strcmp(allNames, 'Series4') | strcmp(allNames, 'Series11')|...
    strcmp(allNames, 'Series5') | strcmp(allNames, 'Series8')|...
    strcmp(allNames, 'Series6') | strcmp(allNames, 'Series9')|...
    strcmp(allNames, 'Series12')| strcmp(allNames, 'Series7');


aB0=   (strcmp(allNames, 'Series13')    )| (strcmp(allNames, 'Series23')    )|...
    strcmp(allNames, 'Series14') | strcmp(allNames, 'Series18')|...
    strcmp(allNames, 'Series17') | strcmp(allNames, 'Series19')|...
    strcmp(allNames, 'Series15') | strcmp(allNames, 'Series21')|...
    strcmp(allNames, 'Series16') | strcmp(allNames, 'Series22')|...
    strcmp(allNames, 'Series20') | strcmp(allNames, 'Series24');



aBNight=   (strcmp(allNames, 'Series25')    )| (strcmp(allNames, 'Series31')    )|...
    strcmp(allNames, 'Series26') | strcmp(allNames, 'Series32')|...
    strcmp(allNames, 'Series27') | strcmp(allNames, 'Series33')|...
    strcmp(allNames, 'Series28') | strcmp(allNames, 'Series34')|...
    strcmp(allNames, 'Series29') | strcmp(allNames, 'Series35')|...
    strcmp(allNames, 'Series30') | strcmp(allNames, 'Series36');


abECD=   (strcmp(allNames, 'Series37')    )| (strcmp(allNames, 'Series43')    )|...
    strcmp(allNames, 'Series38') | strcmp(allNames, 'Series44')|...
    strcmp(allNames, 'Series39') | strcmp(allNames, 'Series45')|...
    strcmp(allNames, 'Series40') | strcmp(allNames, 'Series46')|...
    strcmp(allNames, 'Series41') | strcmp(allNames, 'Series47')|...
    strcmp(allNames, 'Series42') | strcmp(allNames, 'Series48');



randArray=cell(1,length(allNames));

    for tt=1:length(allNames)
randArray{tt}=(allNames{tt}(7:end));
    end

randArray=categorical(randArray);

genoArray=nan(1,length(allNames)); genoArray(isWT)=2; genoArray(isHet)=1; genoArray(isKO)=0; 

treatArray=cell(1,length(allNames)); treatArray(isPBS)={'PBS'}; 
treatArray(aB0)={'aBTP0'}; treatArray(aBNight)={'aBNightB4'};  
treatArray(abECD)={'aB + ECD'};

treatArray=categorical(treatArray);


%% statistical testing

val2Test=avgBodipFluroCellProx;



% Build table
glmeTable = table( val2Test',...
    treatArray',...
    (genoArray)',...
    randArray',...
    'VariableNames', {'testVal','Treatment','GenoType','RandFactor'});
% 
% 
% 
 fitglme(glmeTable, 'testVal~ GenoType*Treatment + (1|RandFactor)')
% 

        % avgBodipFluroCellProx(tt)=nanmean(bodipFluoroCollectNuc); %likewise for aB fluro
        % avgBodipAreaCellProx(tt)=nanmean(bodipAreaCollectNuc); %likewise for aB area
        % 
        % numBodipCellProx(tt)= sum(bodipTrueFalseNuc)/  (length(bodipTrueFalseNuc)); % this will be the ratio of cells with at least one aB object within their radius 
        % 
        % avgBodipFluoroTot(tt)=nanmean(bodipFluoroCollectTot(bodipFluoroCollectTot>= aBFluoroCutoff )); % this will be the number of aB spots in an image (normalized to nuclei count)
        % avgBodipAreaTot(tt)=nanmean(bodipAreaCollectTot(bodipFluoroCollectTot>= aBFluoroCutoff ) ); % this will be the average fluoresence across the whole image
        % numBodipSpots(tt)=  (sum(bodipFluoroCollectTot>= aBFluoroCutoff )   / height(nucFeaturesAtPlay)  ); % likewise with area
%% making a nice histo




% allGeno=[isWT&abECD; isHet&abECD; isKO&abECD];
% figure
% hold on
% for genoCount=1:3
%     bodipyFluoroCollect=[];
%     for tt=  1: find(allGeno(genoCount,:))
%                         rawBodipIm=rawBodipCollect{tt};
%                         bodipFeaturesAtPlay=bodipFeatureCollect{tt};
% 
%                         for zz=1:height(bodipFeaturesAtPlay)
%                             bodipyFluoroCollect(1+length(bodipyFluoroCollect))= nanmean(rawBodipIm(bodipFeaturesAtPlay.PixelIdxList{zz})  ); %simply grab the pixel values found in the space of our mask
% 
%                         end
% 
%     end
% histogram(bodipyFluoroCollect)
% 
% end
% 
% 
% 
% 
% 


%% analysis

val2Test=avgBodipFluroCellProx;
val2Test(isinf(val2Test))=nan;

% val2Test(val2Test >    (nanmean(val2Test)  +4*nanstd(val2Test)  )     )=nan;



ylablr=' Average Bodipy Fluoresence';

titlr=' Average Bodipy Fluoresence PBS vs 24 Hour aB +ECD  ';


treat1=aBNight; treat1Name='aB 24 hours';
treat2=abECD; treat2Name='aB 24 Hours Prior + ECD';




figure
hold on

b=bar(1,nanmean(val2Test(isWT & treat1)))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [.8 .2 .5]; 
ylabel(ylablr,'FontSize',20)
a=gca; a.XTickLabel=[];



b=bar(3,nanmean(val2Test(isHet & treat1)))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0 0.7 .25];


b=bar(5,nanmean(val2Test(isKO & treat1 )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0.3 0.1 .6];


b=bar(7,nanmean(val2Test((isWT & treat2) )   ))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [.9 .1 .1]; 


b=bar(9,nanmean(val2Test(isHet & treat2)))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0.2 0.6 .3];



b=bar(11,nanmean(val2Test(isKO & treat2 )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0.4 0.5 .7];


hold on
a=scatter(rand(1, sum((isWT & treat1)))+.5, val2Test(isWT & treat1), 'Marker', 'o' );
a.CData=[.7 .3 .6]; 

b=scatter(rand(1, sum(isHet & treat1))+2.5, val2Test(isHet & treat1), 'Marker', 'o' );
b.CData(1,:) = [0 0.7 .25];

c=scatter(rand(1, sum(isKO & treat1))+4.5, val2Test(isKO & treat1), 'Marker', 'o' );
c.CData(1,:) = [0.3 0.1 .6];

a=scatter(rand(1, sum((isWT & treat2)))+6.5, val2Test(isWT & treat2), 'Marker', 'o' );
a.CData=[.9 .1 .1]; 


b=scatter(rand(1, sum(isHet & treat2))+8.5, val2Test(isHet & treat2), 'Marker', 'o' );
b.CData(1,:) = [0 0.9 .1];


c=scatter(rand(1, sum(isKO & treat2))+10.5, val2Test(isKO & treat2), 'Marker', 'o' );
c.CData(1,:) = [0.3 0.2 .7];


legend({['WT ',treat1Name],...
    ['HET ',treat1Name],...
    ['KO ',treat1Name],...
    ['WT ',treat2Name],...
    ['HET ',treat2Name],...
    ['KO' ,treat2Name]})



title(titlr,'FontSize',30)



%anovan(val2Test(uM1), {genoArray(uM1)}  )





 % [p, tbl, stats] = kruskalwallis(val2Test(uM1), genoArray(uM1));