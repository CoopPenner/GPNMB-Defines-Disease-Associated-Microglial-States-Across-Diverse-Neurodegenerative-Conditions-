%abeta joint analysis




  set1Data=readtable('/Volumes/PC60/phagocytosisAssays/ABeta_Uptake/allResults_aBetaUptakeSet3/allDataResultsSet1New.xlsx');
  set2Data=readtable('/Volumes/PC60/phagocytosisAssays/ABeta_Uptake/allResults_aBetaUptakeSet3/allDataResultsSet2.xlsx');
    
set2Data.Properties.VariableNames{'avg_aB_Area_Total_Image_normalized'} = 'avg_aB_Area_Total_Image_Normalized';
set2Data.Properties.VariableNames{'avg_aB_Fluoro_Total_Image_normalized'}='avg_aB_Fluoro_Total_Image_Normalized'
  

  set1AllNames=set1Data.ImageCode;
   
    WT= contains(set1AllNames,'_01_')| contains(set1AllNames,'_02_') | contains(set1AllNames,'_03_') ;
    HET=contains(set1AllNames,'_04_')| contains(set1AllNames,'_05_') | contains(set1AllNames,'_06_') ;
    KO1=contains(set1AllNames,'_07_')| contains(set1AllNames,'_08_') | contains(set1AllNames,'_09_') ;
    KO2=contains(set1AllNames,'_10_')| contains(set1AllNames,'_11_') | contains(set1AllNames,'_12_') ;
    

    WT(  contains(set1AllNames,'_01_') &  contains(set1AllNames,'A')   )=false;
    KO2(  contains(set1AllNames,'_10_') &  contains(set1AllNames,'B')   )   =false;
    KO2(  contains(set1AllNames,'_11_') &  contains(set1AllNames,'C')   )   =false;
    PBS= contains(set1AllNames,'A');
    nM200=contains(set1AllNames,'B');
    uM1=contains(set1AllNames,'C');


   
   genoArraySet1=cell(1,length(set1AllNames)); genoArraySet1(WT)={'WT'}; genoArraySet1(HET)={'HET'}; genoArraySet1(KO1)={'KO'}; 
       % Replace empty cells with 'unknown' (or whatever label fits)
        genoArraySet1(cellfun(@isempty, genoArraySet1)) = {'filtered'};
        % Now convert safely
       % genoArraySet1 = categorical(genoArraySet1);
   
  % genoArraySet1=categorical(genoArraySet1);

   treatArraySet1=cell(1,length(set1AllNames)); treatArraySet1(PBS)={'PBS'}; treatArraySet1(nM200)={'nM200'}; treatArraySet1(uM1)={'uM1'}; %treatArraySet1=categorical(treatArraySet1);

    set1Data.Genotype=genoArraySet1';
    set1Data.Treatment=treatArraySet1';

set2DataAllNames=set2Data.ImageCode;

 
     WT=contains(set2DataAllNames,'_10_')| contains(set2DataAllNames,'_11_') | contains(set2DataAllNames,'_12_') ;
     HET=contains(set2DataAllNames,'_04_')| contains(set2DataAllNames,'_05_') | contains(set2DataAllNames,'_06_') ;
     KO1=contains(set2DataAllNames,'_07_')| contains(set2DataAllNames,'_08_') | contains(set2DataAllNames,'_09_') ;
    
     WT(  contains(set2DataAllNames,'_11_') &  contains(set2DataAllNames,'C')   )=false;
    WT(  contains(set2DataAllNames,'_10_') &  contains(set2DataAllNames,'B')   )=false;
    WT(  contains(set2DataAllNames,'_12_') &  contains(set2DataAllNames,'A')   )   =false;
    KO1(  contains(set2DataAllNames,'_01_') &  contains(set2DataAllNames,'A')   )   =false;

 
    PBS= contains(set2DataAllNames,'A');
    nM200=contains(set2DataAllNames,'B');
    uM1=contains(set2DataAllNames,'C');

   
   genoArraySet2=cell(1,length(set2DataAllNames)); genoArraySet2(WT)={'WT'}; genoArraySet2(HET)={'HET'}; genoArraySet2(KO1)={'KO'}; 
       % Replace empty cells with 'unknown' (or whatever label fits)
        genoArraySet2(cellfun(@isempty, genoArraySet2)) = {'filtered'};
        % Now convert safely
        % genoArraySet2 = categorical(genoArraySet2);

   % genoArraySet2=categorical(genoArraySet2);

   treatArraySet2=cell(1,length(set2DataAllNames)); treatArraySet2(PBS)={'PBS'}; treatArraySet2(nM200)={'nM200'}; treatArraySet2(uM1)={'uM1'}; %treatArraySet2=categorical(treatArraySet2);
        
        set2Data.Genotype=genoArraySet2';
        set2Data.Treatment=treatArraySet2';

%%

set1Treat=set1Data.Treatment;

val2TestSet1 = set1Data.Ratio_of_Cells_Pos_for_aB_normalized; set1Geno = set1Data.Genotype;
val2TestSet1(val2TestSet1> nanmean(val2TestSet1)+4*nanstd(val2TestSet1))=nan;


normVal2TestSet1 = val2TestSet1 / nanmean(val2TestSet1(strcmp(set1Geno,'WT') &    strcmp(set1Treat,'nM200')    )   )  ;%(strcmp(set1Treat,'uM1') |


set2Treat=set2Data.Treatment;
val2TestSet2 = set2Data.Ratio_of_Cells_Pos_for_aB_normalized; set2Geno = set2Data.Genotype;
val2TestSet2(val2TestSet2> nanmean(val2TestSet2)+4*nanstd(val2TestSet2))=nan;

normVal2TestSet2 = val2TestSet2 / nanmean(val2TestSet2(strcmp(set2Geno,'WT')& (strcmp(set2Treat,'nM200')       )   ))  ;%| strcmp(set2Treat,'nM200') 

ylablr='Normalized aB Fluorescence (au)';
titlr='Normalized aB Fluorescence Across Genotypes';



set1Data.normVal = normVal2TestSet1;
set2Data.normVal = normVal2TestSet2;

% Add Set label
set1Data.Set = repmat("Set1", height(set1Data), 1);
set2Data.Set = repmat("Set2", height(set2Data), 1);


totTable = [set1Data; set2Data];
totTable.Set = categorical(totTable.Set);

geno = totTable.Genotype; 
treatr=totTable.Treatment;

val2Test = totTable.normVal;


isSet1 = totTable.Set == "Set1";
isSet2 = totTable.Set == "Set2";

WT  = contains(geno,'WT');
ECD = contains(geno,'ECD');
HET = contains(geno,'HET');
KO  = ~contains(geno,'ECD') & ~HET & contains(geno,'KO');


uM1=contains(treatr,'uM1');
nM200=contains(treatr,'nM200');

treatAtPlay=nM200;

% Colors per genotype
wtCol  = [.8 .2 .5];
hetCol = [0 0.7 .25];
koCol  = [0.3 0.1 .6];

filterCatch=strcmp(geno,'filtered');
totTable(filterCatch,:)=[];
    



% Colors per genotype
wtCol  = [.8 .2 .5];
hetCol = [0 0.7 .25];
koCol  = [0.3 0.1 .6];


figure; hold on

% Bars
b = bar(1, nanmean(val2Test(WT&treatAtPlay)));  b.FaceColor='flat'; b.FaceAlpha=.3; b.BarWidth=1.5; b.CData(1,:)=wtCol;
b = bar(3, nanmean(val2Test(HET&treatAtPlay))); b.FaceColor='flat'; b.FaceAlpha=.3; b.BarWidth=1.5; b.CData(1,:)=hetCol;
b = bar(5, nanmean(val2Test(KO&treatAtPlay)));  b.FaceColor='flat'; b.FaceAlpha=.3; b.BarWidth=1.5; b.CData(1,:)=koCol;


val2Test(val2Test> (nanmean(val2Test)+2*nanstd(val2Test)))=nan;

ylabel(ylablr, 'FontSize', 30);
a = gca; a.XTickLabel = [];

title(titlr,'FontSize',30)

% plottin 
groups    = {WT&treatAtPlay , HET&treatAtPlay , KO&treatAtPlay };
xPos      = [.5, 2.5, 4.5, 6.5];
colors    = {wtCol, hetCol, koCol};

for gg = 1:3
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



totTable.Genotype = reordercats(totTable.Genotype, ["WT","HET","KO"]);

totTable.randArray=categorical(totTable.ImageCode);



glme = fitglme(totTable, ...
    'Ratio_of_Cells_Pos_for_aB_normalized ~ 1 + Genotype*Treatment + Set + (1|randArray)', ...
    'Distribution', 'Normal', ...
    'Link', 'identity');
disp(glme)


