%DQBSA Analysis Comp Main


%% load both genotypes
WTData=readtable('/Volumes/PC60/iMicroVids/lysoFunction/7_17_26/WT_DQBSA/WT_DQBSA_bodipyPerFrame_7_17_26.csv');
KOData=readtable('/Volumes/PC60/iMicroVids/lysoFunction/7_17_26/KO_DQBSA/KO_DQBSA_bodipyPerFrame_7_17_26.csv');
WTData.video=string(WTData.video); WTData.genotype=string(WTData.genotype);
KOData.video=string(KOData.video); KOData.genotype=string(KOData.genotype);
bothData=[WTData; KOData];
set1=true;
%% summarize, one row per video per timepoint
vidList=unique(bothData.video);
meanTable=table();
for vv=1:length(vidList)
thisVid=bothData(bothData.video==vidList(vv), :);
frameList=unique(thisVid.frame);
for ff=1:length(frameList)
        thisFrame=thisVid(thisVid.frame==frameList(ff), :);
        nObj=height(thisFrame);
        meanTable=[meanTable; table( vidList(vv), ...
            unique(thisFrame.genotype), ...
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
'VariableNames', {'video','genotype','frame','timeMin','nObj','meanFluoro','stdFluoro','semFluoro','meanArea','semArea','totArea','meanIntDen','totIntDen','meanFluoroMax'})];
end
end


%% drop spurious KO video 4 timepoints (debris occlusion at 75 and 90 min)in only set 1

if set1

badRows = meanTable.genotype=="KO" & ...
          contains(meanTable.video, "_M4") & ...
          ismember(meanTable.timeMin, [75 90]);
meanTable(badRows, :) = [];

end


%% plot, KO in greens, WT in pinks, with SEM bands
koVids=unique(meanTable.video(meanTable.genotype=="KO"));
wtVids=unique(meanTable.video(meanTable.genotype=="WT"));
% shade ramps, dark to light within each genotype
greenMap=[linspace(0.00,0.55,length(koVids))', linspace(0.35,0.80,length(koVids))', linspace(0.10,0.45,length(koVids))'];
pinkMap =[linspace(0.55,0.95,length(wtVids))', linspace(0.05,0.55,length(wtVids))', linspace(0.35,0.70,length(wtVids))'];
figure; hold on
for vv=1:length(koVids)
    thisVid=meanTable(meanTable.video==koVids(vv), :);
    xFill=[thisVid.timeMin; flipud(thisVid.timeMin)];
    yFill=[thisVid.meanFluoro+thisVid.semFluoro; flipud(thisVid.meanFluoro-thisVid.semFluoro)];
    fill(xFill, yFill, greenMap(vv,:), 'FaceAlpha', .15, 'EdgeColor','none', 'HandleVisibility','off');
    plot(thisVid.timeMin, thisVid.meanFluoro, '-o', 'Color', greenMap(vv,:), ...
'LineWidth', 1.5, 'MarkerSize', 4, 'HandleVisibility','off');
end
for vv=1:length(wtVids)
    thisVid=meanTable(meanTable.video==wtVids(vv), :);
    xFill=[thisVid.timeMin; flipud(thisVid.timeMin)];
    yFill=[thisVid.meanFluoro+thisVid.semFluoro; flipud(thisVid.meanFluoro-thisVid.semFluoro)];
    fill(xFill, yFill, pinkMap(vv,:), 'FaceAlpha', .15, 'EdgeColor','none', 'HandleVisibility','off');
    plot(thisVid.timeMin, thisVid.meanFluoro, '-o', 'Color', pinkMap(vv,:), ...
'LineWidth', 1.5, 'MarkerSize', 4, 'HandleVisibility','off');
end
% dummy handles so the legend shows only genotype
hKO=plot(nan, nan, '-o', 'Color', greenMap(round(end/2),:), 'LineWidth', 2, 'DisplayName','KO');
hWT=plot(nan, nan, '-o', 'Color', pinkMap(round(end/2),:), 'LineWidth', 2, 'DisplayName','WT');
xlabel('time (min)'); ylabel('mean DQ-BSA fluoresence per object')
title('DQ-BSA over time')
legend([hWT hKO], 'Location','eastoutside'); box off
%% glme
meanTable.genotype=categorical(meanTable.genotype);
meanTable.genotype=reordercats(meanTable.genotype, {'WT','KO'}); % WT as reference
meanTable.video=categorical(meanTable.video);
glmeMod=fitglme(meanTable, 'meanFluoro ~ 1 + genotype*timeMin + nObj + (1|video)');
disp(glmeMod)
disp(anova(glmeMod))