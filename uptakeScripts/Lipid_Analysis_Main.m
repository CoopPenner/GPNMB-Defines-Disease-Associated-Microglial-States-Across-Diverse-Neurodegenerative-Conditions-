%basic image analysis for Bodipy analysis

%% 1) loading everything and setting initial parameters
%load images (these are spearated into three channels, C1 Trem2 C2 aB C3
%Nuclei)
%this is designed for a 20x field of view. 

addpath('/Users/pennerc/Documents/Documents - BGS-STU23-138/MATLAB/GPNMB_Global_Analysis/GPNMB-Global-Analysis');

set1=false;
set2=false;
set3=false;
Well24=true;
%Link the folder containing all image 

if set1
homePath='/Volumes/PC60/phagocytosisAssays/Lipids/Set1_96/RowG_BODIPY-FL';
elseif set2
  homePath=  '/Volumes/PC60/phagocytosisAssays/Lipids/Set2/Bodipy_FL';
elseif set3
   homePath='/Volumes/PC60/phagocytosisAssays/Lipids/Set3_96/BODIPY_FL';

elseif Well24
 homePath='/Volumes/PC60/phagocytosisAssays/Lipids/24Well/untouched tifs';
end
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
    lipidCircCut=0; % this is a cutoff for the minimum circularity of a detected object, less circular objects tend to be debris n whatnot.  
    lipidMaxSize=30000; % 
    lipidMinSize=25;
    lipidStrelType='arbitrary'; %don't want to bias shape of aB 
    lipidStrelRadiusSize=3;



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
if set1
   lipidFluoroCutoff=10000; %this is the minimum fluoro for an object to be considered "real" It was calaculated through visual assesment of the fluoro distribution
   % this was output in section 3 below, there is a trimodal distribution
   % with a large peak of weird noise at around 7000

elseif set2 | set3 | Well24
   lipidFluoroCutoff=100; 
   lipidFluoroCutoff=100; 

end
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
lipidMaskCollect=cell(1,length(allNames)); lipidFeatureCollect= cell(1,length(allNames)); rawLipidCollect= cell(1,length(allNames));
T2MaskCollect=cell(1,length(allNames)); T2FeatureCollect= cell(1,length(allNames)); rawT2Collect= cell(1,length(allNames));




%% 2) capturing basic masks for all channels as well as raw data
    %iterating through each name. this will contain three different
    %channels this will take a while



%initializing wait bar
f = waitbar(0, 'Starting');

  
 for dd= 1:length(allNames)
         %I'm fairly certain I could just iterate through in sets of 3, but
         %I want to be agnostic to ordering in case that's not always true.
     nameBaseAtPlay=allNames{dd};

     if set1
      nuclearImAtPlay=imread( [homePath,'/',nameBaseAtPlay,'_DRAQ5.tif']); rawNucCollect{dd}=nuclearImAtPlay; % I grab the raw images for downsttream plotting and intensity filtering
      Trem2ImAtPlay=imread( [homePath,'/',nameBaseAtPlay,'_TREM2.tif']);   rawT2Collect{dd}=Trem2ImAtPlay;

     lipidImAtPlay=imread( [homePath,'/',nameBaseAtPlay,'_Lipid.tif']);      rawLipidCollect{dd}=lipidImAtPlay;

     elseif set2
            nuclearImAtPlay=imread( [homePath,'/',nameBaseAtPlay,'_DRAQ5.tif']); rawNucCollect{dd}=nuclearImAtPlay; % I grab the raw images for downsttream plotting and intensity filtering
             Trem2ImAtPlay=imread( [homePath,'/',nameBaseAtPlay,'_TREM2.tif']);   rawT2Collect{dd}=Trem2ImAtPlay;
             lipidImAtPlay=imread( [homePath,'/',nameBaseAtPlay,'_BODIPY.tif']);      rawLipidCollect{dd}=lipidImAtPlay;
     elseif set3

            nuclearImAtPlay=imread( [homePath,'/',nameBaseAtPlay,'_DRAQ5.tif']); rawNucCollect{dd}=nuclearImAtPlay; % I grab the raw images for downsttream plotting and intensity filtering
             Trem2ImAtPlay=imread( [homePath,'/',nameBaseAtPlay,'_TREM.tif']);   rawT2Collect{dd}=Trem2ImAtPlay;
             lipidImAtPlay=imread( [homePath,'/',nameBaseAtPlay,'_BODIPY-FL.tif']);      rawLipidCollect{dd}=lipidImAtPlay;


     
     elseif Well24
            nuclearImAtPlay=imread( [homePath,'/',nameBaseAtPlay,'_DAPI.tif']); rawNucCollect{dd}=nuclearImAtPlay; % I grab the raw images for downsttream plotting and intensity filtering
             Trem2ImAtPlay=nan;   rawT2Collect{dd}=nan;
        lipidImAtPlay=imread( [homePath,'/',nameBaseAtPlay,'_BODIPY.tif']);      rawLipidCollect{dd}=lipidImAtPlay;



     end

    [nucMaskCollect{dd},~,nucFeatureCollect{dd}] = aBNuclearMaskMaker(nuclearImAtPlay,nucMinSize, nucMaxSize, nuccircCut,nucStrelRadiusSize,nucStrelType,demoPlot);

    
if demoPlot
  title('nuclear')  
end





if set2
    c=true;
    backgroundGaussian=60;
end

    [lipidMaskCollect{dd},~,lipidFeatureCollect{dd}] = lipidMaskMaker(lipidImAtPlay,lipidMinSize, lipidMaxSize, lipidCircCut,lipidStrelRadiusSize,lipidStrelType,backgroundCorrect,backgroundGaussian,demoPlot);


if demoPlot
      title('Lipid')  

end

if ~Well24
            [T2MaskCollect{dd},~,T2FeatureCollect{dd}] = aBNuclearMaskMakerTrem2Specific(Trem2ImAtPlay,T2minSize, T2maxSize, T2CircCut,T2strelRadiusSize,T2strelType, T2detectSensitivity, demoPlot);
        
        if demoPlot
        
            title('trem2')  
            
        end
end
    waitbar(dd/length(allNames), f, sprintf('Progress: %d %%', floor( (dd/length(allNames))  *100)));
    pause(0.1);



     totIm = cat(3,  false(size(lipidImAtPlay)),lipidImAtPlay, nuclearImAtPlay);     % 3D array HxWx3
     maskImTot = cat(3,false(size(lipidImAtPlay)), lipidMaskCollect{dd}, nucMaskCollect{dd});  
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



%% 3) loading all Fluoresences to generate a cutoff (you can do this automatically with just some kind of std metric but I think it's better to finesse it)


lipidFluoroCollect=[];
if set1 || set2
KO=contains(allNames,'KO') & ~contains(allNames,'ECD');
WT=contains(allNames,'WT');
ECD=contains(allNames,'ECD');

HET=contains(allNames,'HET');

idx=find(ECD);

elseif Well24


WT=contains(allNames,'Series1') | contains(allNames,'Series2') | contains(allNames,'Series3') |...
    contains(allNames,'Series4') | contains(allNames,'Series13') | contains(allNames,'Series14') |...
    contains(allNames,'Series25') | contains(allNames,'Series26') | contains(allNames,'Series27') |...
    contains(allNames,'Series28') | contains(allNames,'Series37') | contains(allNames,'Series38') |...
    contains(allNames,'Series39') | contains(allNames,'Series40') ;


KO=contains(allNames,'Series9') | contains(allNames,'Series10') | contains(allNames,'Series11') |...
    contains(allNames,'Series12') | contains(allNames,'Series21') | contains(allNames,'Series22') |...
    contains(allNames,'Series23') | contains(allNames,'Series24') | contains(allNames,'Series33') |...
    contains(allNames,'Series34') | contains(allNames,'Series35') | contains(allNames,'Series36') |...
    contains(allNames,'Series45') | contains(allNames,'Series46') | contains(allNames,'Series47') | contains(allNames,'Series48') ;

HET=contains(allNames,'Series7') | contains(allNames,'Series8') | contains(allNames,'Series5') |...
    contains(allNames,'Series6') | contains(allNames,'Series19') | contains(allNames,'Series20') |...
    contains(allNames,'Series29') | contains(allNames,'Series30') | contains(allNames,'Series31') |...
    contains(allNames,'Series32') | contains(allNames,'Series41') | contains(allNames,'Series42') |...
    contains(allNames,'Series43') | contains(allNames,'Series44') | contains(allNames,'Series45') | contains(allNames,'Series46') ;


WT=contains(allNames,'Series1') | contains(allNames,'Series2') | contains(allNames,'Series3') |...
    contains(allNames,'Series4') | contains(allNames,'Series13') | contains(allNames,'Series14') |...
    contains(allNames,'Series25') | contains(allNames,'Series26') | contains(allNames,'Series27') |...
    contains(allNames,'Series28') | contains(allNames,'Series37') | contains(allNames,'Series38') |...
    contains(allNames,'Series39') | contains(allNames,'Series40') ; WT(KO)=false; WT(HET)=false;






end


idx=find(HET);
            for tt=1:length(idx)
                    rawLipidIm= rawLipidCollect{idx(tt)};
                    % rawT2Im=rawT2Collect{tt};
                    lipidFeaturesAtPlay=lipidFeatureCollect{idx(tt)};
                    % T2FeaturesAtPlay=T2FeatureCollect{tt};

                    for zz=1:height(lipidFeaturesAtPlay)
                        lipidFluoroCollect(1+length(lipidFluoroCollect))= nanmean(rawLipidIm(lipidFeaturesAtPlay.PixelIdxList{zz})  ); %simply grab the pixel values found in the space of our mask

                    end


            end

% based on just looking it seems like a Fluoro Cut off of 10 is reasonable
% for both channels. There is a clear high degree of low signal Trem2, less
% so for aB. 

             % 
             % figure
             % histogram(T2FluoroCollect)

             
hold on

histogram(lipidFluoroCollect)

nanmean(lipidFluoroCollect)
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



            % we will iterate through each name, and then anchor analyses
            % on a per nucleus basis. 

        numLipidSpotsTot=nan(1,length(allNames)); % this will be the number of aB spots within our radius of a nuclei
        avgLipidFluoroTot=nan(1,length(allNames)); %likewise for aB fluro
        avgAreaLipidTot=nan(1,length(allNames)); %likewise for aB area
        areaLipidTotCollect=cell(1,length(allNames)); %likewise for aB area
        fluoroLipidTotCollect=cell(1,length(allNames)); %likewise for aB area




       for tt=1: length(allNames)

                nucFeaturesAtPlay=nucFeatureCollect{tt}; % loop through nuclei in a given image
    
                rawLipidIm= rawLipidCollect{tt}; % to grab fluoresence 
                lipidFeaturesAtPlay=lipidFeatureCollect{tt};% for size and location
                lipidMaskAtPlay=lipidMaskCollect{tt}; % for mapping
                lipidDistMap = bwdist(lipidMaskAtPlay, 'euclidean');
    

                lipidAreaCollectTot=nan(1,height(nucFeaturesAtPlay));
                lipidFluoroCollectTot=nan(1,height(nucFeaturesAtPlay));

                validIdx=false(1,height(lipidFeaturesAtPlay));
                lipidFluoroCollect=nan(1,height(lipidFeaturesAtPlay));

                                %first grab all aB objects within 2 pixels of a given nucleus (anchored to cells) 
                      for dd=1:height(lipidFeaturesAtPlay)


                          idxAtPlay= lipidFeaturesAtPlay.PixelIdxList{dd};

                          lipidy=rawLipidIm(idxAtPlay);
                          lipidFluoroCollect(dd)=nanmean(lipidy);

                          if nanmean(lipidy) >= lipidFluoroCutoff
                              validIdx(dd)=true;
                          else
                              validIdx(dd)=false;
                          end
                   


                      end
                allArea=lipidFeaturesAtPlay.Area(validIdx);
         


        numLipidSpotsTot(tt)=sum(validIdx); % this will be the number of aB spots within our radius of a nuclei
        avgLipidFluoroTot(tt)=nanmean(lipidFluoroCollect(validIdx)); %likewise for aB fluro
        avgAreaLipidTot(tt)=nanmean(allArea); %likewise for aB area

        areaLipidTotCollect{tt}=allArea; %likewise for aB area
        fluoroLipidTotCollect{tt}=lipidFluoroCollect(validIdx); %likewise for aB area



          end 


            

isAbNightBefore = double(abNightBefore(:));
isAddECD        = double(addECD(:));

% Combine into table — keep only WT/HET/KO rows
tbl = table(numLipidSpotsTot(:),...
    avgLipidFluoroTot(:),...
    avgAreaLipidTot(:),...
    areaLipidTotCollect(:),...
    fluoroLipidTotCollect(:),...
    genotype(:),...
    isabNightBefore(:),...
    isAddECD(:),...
    'VariableNames', {'numLipidSpotsTot','avgLipidFluoroTot','avgAreaLipidTot','areaLipidTotCollect','fluoroLipidTotCollect','genotype','aBNightBefore','addECD'});


writetable(tbl,'/Users/pennerc/Desktop/scRNASeqScriptsCKP/24WellQuantsCKP.csv')

%% plotting and stats    

val2Test=avgLipidFluoroTot;
val2Test(isinf(val2Test))=nan;




ylablr=' Average Fluoresence of all Detected Bodipy-C16 Objects (au)';

titlr=' BODIPY-C16 Fluoresence  ';







if set1  

 WT=contains(allNames,'WT');
 ECD=contains(allNames,'ECD');
 HET= contains(allNames, 'WellG04')|contains(allNames, 'WellG05')|contains(allNames, 'WellG06');
 KO= ~contains(allNames,'ECD') & ~HET &...
     contains(allNames, 'WellG07')|contains(allNames, 'WellG08')|contains(allNames, 'WellG09');



elseif set2


 WT=contains(allNames,'WT');
 ECD=contains(allNames,'ECD');
 HET= contains(allNames, 'HET');
 KO= ~contains(allNames,'ECD') & ~HET &...
     contains(allNames, 'KO');

normalizedVal2Test=val2Test/ nanmean(val2Test(WT));

figure
hold on

b=bar(1,nanmean(val2Test(WT)))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [.8 .2 .5]; 
ylabel(ylablr,'FontSize',30)
a=gca; a.XTickLabel=[];



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


b=bar(7,nanmean(val2Test(ECD )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0 0.2 .9];


hold on
a=scatter(rand(1, sum((WT & ~ECD )))+.5, val2Test(WT & ~ECD), "filled", 'Marker', 'o' );
a.CData=[.7 .2 .4]; 

b=scatter(rand(1, sum(HET & ~ECD ))+2.5, val2Test(HET & ~ECD ), "filled", 'Marker', 'o'  );
b.CData(1,:) = [0 0.7 .25];


c=scatter(rand(1, sum(KO& ~ECD ) )+4.5, val2Test(KO& ~ECD ), "filled", 'Marker', 'o' );
c.CData(1,:) = [0.3 0.1 .6];

d=scatter(rand(1, sum(ECD ) )+6.5, val2Test(ECD ), "filled", 'Marker', 'o' );
d.CData(1,:) = [0 0.2 .9];






% Genotype labels (categorical)
genotype = repmat("WT", length(val2Test), 1);
genotype(HET) = "HET";
genotype(KO)  = "KO";
genotype(ECD)= "KO+ECD";
genotype = categorical(genotype, ["WT","HET","KO","KO+ECD"]);  % WT as reference

% Binary predictors
isAddECD        = double(ECD(:));

imageNum = categorical(allNames);



tbl = table(val2Test(:), genotype, isAddECD, imageNum', normalizedVal2Test', ...
    'VariableNames', {'BodipyFluoro','Genotype','addECD','ImageNum','NormalizedFluoro' });

glme = fitglme(tbl, ...
    'BodipyFluoro ~ 1 + Genotype + (1|ImageNum)', ...
    'Distribution', 'Normal', ...   % change if val2Test is not continuous/normal
    'Link', 'identity');            % change accordingly (e.g. 'log' for Poisson)

disp(glme)
if set1
writetable(tbl, ['/Volumes/PC60/lipidStuff','/','set1Results.xls' ])
elseif set2
  writetable(tbl, ['/Volumes/PC60/lipidStuff','/','set2Results.xls' ])
end

elseif Well24
WT=contains(allNames,'Series1') | contains(allNames,'Series2') | contains(allNames,'Series3') |...
    contains(allNames,'Series4') | contains(allNames,'Series13') | contains(allNames,'Series14') |...
    contains(allNames,'Series25') | contains(allNames,'Series26') | contains(allNames,'Series27') |...
    contains(allNames,'Series28') | contains(allNames,'Series37') | contains(allNames,'Series38') |...
    contains(allNames,'Series39') | contains(allNames,'Series40') ;


KO=contains(allNames,'Series9') | contains(allNames,'Series10') | contains(allNames,'Series11') |...
    contains(allNames,'Series12') | contains(allNames,'Series21') | contains(allNames,'Series22') |...
    contains(allNames,'Series23') | contains(allNames,'Series24') | contains(allNames,'Series33') |...
    contains(allNames,'Series34') | contains(allNames,'Series35') | contains(allNames,'Series36') |...
    contains(allNames,'Series45') | contains(allNames,'Series46') | contains(allNames,'Series47') | contains(allNames,'Series48') ;

HET=contains(allNames,'Series7') | contains(allNames,'Series8') | contains(allNames,'Series5') |...
    contains(allNames,'Series6') | contains(allNames,'Series19') | contains(allNames,'Series20') |...
    contains(allNames,'Series29') | contains(allNames,'Series30') | contains(allNames,'Series31') |...
    contains(allNames,'Series32') | contains(allNames,'Series41') | contains(allNames,'Series42') |...
    contains(allNames,'Series43') | contains(allNames,'Series44')  ;
WT(KO)=false; WT(HET)=false;




 abNightBefore=contains(allNames,'Series25') | contains(allNames,'Series26') | contains(allNames,'Series27') |...
     contains(allNames,'Series28') | contains(allNames,'Series29') | contains(allNames,'Series30') |...
     contains(allNames,'Series31') | contains(allNames,'Series32') | contains(allNames,'Series33') |...
     contains(allNames,'Series34') | contains(allNames,'Series35') | contains(allNames,'Series36') ;

addECD=contains(allNames,'Series37') | contains(allNames,'Series38') | contains(allNames,'Series39') |...
    contains(allNames,'Series40') | contains(allNames,'Series41') | contains(allNames,'Series42') |...
    contains(allNames,'Series43') | contains(allNames,'Series44') | contains(allNames,'Series45') |...
    contains(allNames,'Series46') | contains(allNames,'Series47') | contains(allNames,'Series48') ;


normalizedVal2Test=val2Test/ nanmean(val2Test(WT));



figure
hold on

b=bar(1,nanmean(val2Test(WT & ~addECD & ~abNightBefore)))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [.8 .2 .5]; 
ylabel(ylablr,'FontSize',30)
a=gca; a.XTickLabel=[];


b=bar(3,nanmean(val2Test(WT  & addECD & ~abNightBefore)))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [.8 .4 .6]; 
ylabel(ylablr,'FontSize',30)
a=gca; a.XTickLabel=[];




b=bar(5,nanmean(val2Test(HET & ~addECD& ~abNightBefore )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0 0.7 .25];


b=bar(7,nanmean(val2Test(HET & addECD& ~abNightBefore )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0.1 0.7 .4];


b=bar(9,nanmean(val2Test(KO & ~addECD& ~abNightBefore )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0.3 0.1 .6];

b=bar(11,nanmean(val2Test(KO & addECD& ~abNightBefore )))  ;
b.FaceColor = 'flat';
b.FaceAlpha=.3;
b.BarWidth=1.5;
b.CData(1,:) = [0.3 0.1 .6];




hold on
a=scatter(rand(1, sum((WT & ~addECD& ~abNightBefore )))+.5, val2Test(WT & ~addECD& ~abNightBefore  ), "filled", 'Marker', 'o' );
a.CData=[.7 .2 .4]; 

aa=scatter(rand(1, sum((WT & addECD& ~abNightBefore )))+2.5, val2Test(WT & addECD& ~abNightBefore  ), "filled", 'Marker', '<'  );
aa.CData=[.7 .2 .4]; 

% aaa=scatter(rand(1, sum((WT & ~addECD & abNightBefore)))+.5, val2Test(WT & ~addECD & abNightBefore ), "filled", 'Marker',  "square" );
% aaa.CData=[.7 .2 .4]; 

b=scatter(rand(1, sum(HET & ~addECD & ~abNightBefore ))+4.5, val2Test(HET & ~addECD& ~abNightBefore  ), "filled", 'Marker', 'o'  );
b.CData(1,:) = [0 0.7 .25];

bb=scatter(rand(1, sum((HET & addECD& ~abNightBefore )))+6.5, val2Test(HET & addECD& ~abNightBefore  ), "filled", 'Marker', '<' );
bb.CData(1,:) = [0 0.7 .25];

% bbb=scatter(rand(1, sum((HET & ~addECD & abNightBefore)))+2.5, val2Test(HET & ~addECD & abNightBefore ), "filled", 'Marker',  "square" );
% bbb.CData(1,:) = [0 0.7 .25];


c=scatter(rand(1, sum(KO & ~addECD& ~abNightBefore ) )+8.5, val2Test(KO & ~addECD& ~abNightBefore ), "filled", 'Marker', 'o' );
c.CData(1,:) = [0.3 0.1 .6];


cc=scatter(rand(1, sum((KO & addECD& ~abNightBefore )))+10.5, val2Test(KO & addECD& ~abNightBefore ), "filled", 'Marker', '<'  );
cc.CData(1,:) = [0.3 0.1 .6];

% ccc=scatter(rand(1, sum((KO & ~addECD & abNightBefore)))+4.5, val2Test(KO & ~addECD & abNightBefore ), "filled", 'Marker', "square"  );
% ccc.CData(1,:) = [0.3 0.1 .6];

% a=scatter(rand(1, sum((WT )))+6.5, val2Test(WT ), 'Marker', 'o' );
% a.CData=[.9 .1 .1]; 


%Build the data table

% Genotype labels (categorical)
genotype = repmat("WT", length(val2Test), 1);
genotype(HET & ~addECD) = "HET";
genotype(KO & ~addECD)  = "KO";
genotype(WT & ~addECD)  = "WT";
genotype(HET & addECD) = "HET+ECD";
genotype(KO & addECD)  = "KO+ECD";
genotype(WT & addECD)  = "WT+ECD";

genotype = categorical(genotype, ["WT","HET","KO","WT+ECD","HET+ECD","KO+ECD"]);  % WT as reference

% Binary predictors
isAbNightBefore = double(abNightBefore(:));
isAddECD        = double(addECD(:));

imageNum = categorical(allNames);

% Combine into table — keep only WT/HET/KO rows
tbl = table(val2Test(:), genotype, isAddECD, imageNum', ...
    'VariableNames', {'val2Test','Genotype','addECD','ImageNum'});

% validRows = tbl.Genotype == "WT" | tbl.Genotype == "HET" | tbl.Genotype == "KO";
% tbl = tbl(validRows, :);


% Fit the GLME
% Fixed effects: genotype, abNightBefore, addECD (and their interactions if desired)
% Random effect: intercept per image number

glme = fitglme(tbl, ...
    'val2Test ~ 1 + Genotype  + (1|ImageNum)', ...
    'Distribution', 'Normal', ...   % change if val2Test is not continuous/normal
    'Link', 'identity');            % change accordingly (e.g. 'log' for Poisson)

disp(glme)










end




% legend({'WT','HET','KO'},'FontSize',30)
  



title(titlr,'FontSize',30)

%% reload and analyze everything

set1Data = readtable('/Volumes/PC60/phagocytosisAssays/Lipids/set1Results.xls');
set2Data = readtable('/Volumes/PC60/phagocytosisAssays/Lipids/set2Results.xls');

fluorSet1 = set1Data.val2Test; set1Geno = set1Data.Genotype;
normFluorSet1 = fluorSet1 / nanmean(fluorSet1(strcmp(set1Geno,'WT')));
fluorSet2 = set2Data.val2Test; set2Geno = set2Data.Genotype;
normFluorSet2 = fluorSet2 / nanmean(fluorSet2(strcmp(set2Geno,'WT')));

set1Data.normVal = normFluorSet1;
set2Data.normVal = normFluorSet2;

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

% plottin 
groups    = {WT & ~ECD, HET & ~ECD, KO & ~ECD, ECD};
xPos      = [.5, 2.5, 4.5, 6.5];
colors    = {wtCol, hetCol, koCol, ecdCol};

for g = 1:4
    mask = groups{g};
    col  = colors{g};

    % Set 1 
    idx1 = mask & isSet1;
    scatter(rand(1, sum(idx1)) + xPos(g), val2Test(idx1), ...
        "filled", 'Marker', 'o', 'CData', repmat(col, sum(idx1), 1));

    % Set 2 
    idx2 = mask & isSet2;
    scatter(rand(1, sum(idx2)) + xPos(g), val2Test(idx2), ...
        "filled", 'Marker', 'square', 'CData', repmat(col, sum(idx2), 1));
end

% GLME 
totTable.Genotype = categorical(totTable.Genotype, ["WT","HET","KO","KO+ECD"]);

glme = fitglme(totTable, ...
    'normVal ~ 1 + Genotype + Set + (1|ImageNum)', ...
    'Distribution', 'Normal', ...
    'Link', 'identity');
disp(glme)




%% Cemetery


% figure
% 
% valAtPlay=avgAbFluoroCellProx;
% 
% bar([nanmean(valAtPlay(WT & nM200)  ), nanmean(valAtPlay(HET & nM200)  ), nanmean(valAtPlay(  )   ])
% KO & nM200)
% 
% 
% figure
% 
% 
% bar([nanmean(valAtPlay(WT & uM1)  ), nanmean(valAtPlay(HET & uM1)  ), nanmean(valAtPlay(KO & uM1)  )   ])
% 
% 
% 
% 
% figure
% 
% 
% bar([nanmean(valAtPlay(WT & uM200_555)  ), nanmean(valAtPlay(HET & uM200_555)  ), nanmean(valAtPlay(KO & uM200_555)  )   ])



     


            % 
% 
% if testTime

% figure
%  mapperSupreme= blankMap; mapperSupreme(nucFeaturesAtPlay.PixelIdxList{dd})=true;
% maskImTot = cat(3, abMaskAtPlay, mapperSupreme, nucMaskCollect{tt});  
% 
% imshowpair(maskImTot,rawABIm,'montage')
% 
% title(['nucleus ', num2str(dd), ' abIDx=', num2str(sum(abIdx))])
% end

% 
% b=1;
        