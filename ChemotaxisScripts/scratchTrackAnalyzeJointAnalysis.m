%scratchTrack Analyze joint! 

    % 'run_7_6_26',  'WT', '/Volumes/PC60/iMicroVids/setsPassedQC/PBS/WT_KO_PBS_7_6_26/WT/ResultsRepo'
    % 'run_7_6_26',  'KO', '/Volumes/PC60/iMicroVids/setsPassedQC/PBS/WT_KO_PBS_7_6_26/KO/ResultsRepo'

addpath('/Users/pennerc/Documents/Documents - BGS-STU23-138/GitHub/GPNMB_GlobalANalysis_Updated')


genotypes = {'KO', 'WT'};

imFreq = 5;
treatmentPrefix = 'PBS';

% ---- run definitions: {runID, genotype, path} ----
runDefs = {
    'run_6_6_25',  'WT', '/Volumes/PC60/iMicroVids/setsPassedQC/PBS/wt_ko_PBS_6_6_25/WT_PBS/ResultsRepo'
    'run_6_6_25',  'KO', '/Volumes/PC60/iMicroVids/setsPassedQC/PBS/wt_ko_PBS_6_6_25/KO_PBS/ResultsRepo'
    'run_6_21_26', 'WT', '/Volumes/PC60/iMicroVids/setsPassedQC/PBS/WTKOScratchAssayPBSOnly_6_21_26_allChannels/WT/ResultsRepo'
    'run_6_21_26', 'KO', '/Volumes/PC60/iMicroVids/setsPassedQC/PBS/WTKOScratchAssayPBSOnly_6_21_26_allChannels/KO/ResultsRepo'
    'run_7_6_26',  'WT', '/Volumes/PC60/iMicroVids/setsPassedQC/PBS/WT_KO_PBS_7_6_26/WT/ResultsRepo'
     'run_7_6_26',  'KO', '/Volumes/PC60/iMicroVids/setsPassedQC/PBS/WT_KO_PBS_7_6_26/KO/ResultsRepo'

};
%for the time being these are all of our loaded files, I used claude to
%rework my initial scratchTrackAnalyze code, this is a bit less readable.
%But I'll just comment out runDefs between
structAtPlay = struct([]);

for dd = 1:size(runDefs,1)
    runID      = runDefs{dd,1};
    geno       = runDefs{dd,2};
    folderPath = runDefs{dd,3};

    if ~isfolder(folderPath)
        warning('Missing folder: %s', folderPath); continue
    end

    fileList = dir(fullfile(folderPath, '*.mat'));
    fileList = fileList(~startsWith({fileList.name}, '._'));
    fileList = fileList(contains({fileList.name}, 'scratchMetrics'));
    fprintf('%s | %s : %d files\n', runID, geno, numel(fileList));

    for f = 1:numel(fileList)
        loaded = load(fullfile(folderPath, fileList(f).name));
        m = loaded.scratchMet;

        tok = regexp(fileList(f).name, '\((\d+)\)', 'tokens');
        vidNumRaw = str2double(tok{1});

        wellPlate = ceil(vidNumRaw/2);   % 1,2->1  3,4->2  5,6->3  7,8->4 

        structAtPlay(end+1).Genotype   = geno;
        structAtPlay(end).Treatment    = treatmentPrefix;
        structAtPlay(end).runID        = runID;
        structAtPlay(end).VideoNumRaw  = vidNumRaw;
        structAtPlay(end).VideoNum     = sprintf('%s_%s_%d', runID, geno, vidNumRaw);
        structAtPlay(end).wellPlate = sprintf('Well%d_%s', wellPlate, runID);

        structAtPlay(end).InScratchCount   = m.inScratchCount;
        structAtPlay(end).NearScratchCount = m.nearScratchCount;
        structAtPlay(end).TotObjects       = m.totObjects;
        structAtPlay(end).AvgDistToScratch = m.avgDistTot;
        structAtPlay(end).StdDistToScratch = m.stdDistTot;
        structAtPlay(end).ScratchArea      = m.scratchSize;
        structAtPlay(end).wellPlate = sprintf('Well%d_%s_%s', wellPlate, geno, runID);
    end
end



%% building our model data

hoursToPlot = 8;
nFrames     = hoursToPlot*60/imFreq;   % 96 with imFreq=5
ylablr        = 'In Scratch Count';
measureAtPlay = 'scratchCount'; %scratchCount or GlobalDist
rows = {};
for ii = 1:numel(structAtPlay)
    nCells = structAtPlay(ii).TotObjects(:);
    switch measureAtPlay
        case 'globalDist'
            trace = structAtPlay(ii).AvgDistToScratch(:);
        case 'scratchCount'
            trace = structAtPlay(ii).NearScratchCount(:);
            if sum(trace < 0) > 10
                trace = nan(size(trace));
            end
            trace = trace ./ nCells;            % ratio in scratch
    end
    if numel(trace) < nFrames
        warning('%s: only %d frames, skipping', structAtPlay(ii).VideoNum, numel(trace));
        continue
    end
    trace  = trace(1:nFrames);
    nCells = nCells(1:nFrames);
    if ~isfinite(trace(1)) || trace(1) == 0
        warning('%s: bad baseline, skipping', structAtPlay(ii).VideoNum); % I am a little uncertain about these traces that start with zero cells in scratch
        %thus in the current workflow I remove them, they are all in the
        %final replicate and account for 3 total videos (2 in WT and 1 in
        %KO). Removing these does not affect the results. 
        continue
    end
    base     = trace(1);
    distNorm = log(trace ./ base);              % 0 at t=0, log-ratio scale
    distNorm(~isfinite(distNorm)) = NaN;
    t = (0:nFrames-1)' * imFreq;
    rows{end+1} = table( ...
        trace, ...
        distNorm, ...
        repmat(base, nFrames, 1), ...
        nCells, ...
        repmat(structAtPlay(ii).ScratchArea(1), nFrames, 1), ...
        t, ...
        repmat({structAtPlay(ii).Genotype},  nFrames, 1), ...
        repmat({structAtPlay(ii).VideoNum},  nFrames, 1), ...
        repmat({structAtPlay(ii).wellPlate}, nFrames, 1), ...
        repmat({structAtPlay(ii).runID},     nFrames, 1), ...
        'VariableNames', {'Dist','DistNorm','startVal','TotCells','scratchArea','time', ...
                          'genotype','vidNum','wellPlate','runID'});
end
modelData = vertcat(rows{:});
modelData.genotype  = reordercats(categorical(modelData.genotype), {'KO','WT'});
modelData.vidNum    = categorical(modelData.vidNum);
modelData.wellPlate = categorical(modelData.wellPlate);
modelData.runID     = categorical(modelData.runID);
modelData.DistNormClamp = max(modelData.DistNorm, 0);   % negatives -> 0, NaN stays NaN
summary(modelData)



% for the scratch Count we will use the clamped values, ostensibly, in some
% cases when no migration into the scratch occurs the value can come out
% zero, I was worried this would artificially inflate our value... it comes
% out signiicant either way but I think this is more honest. 

%For the global distance we don't use a clamped model because ther can
%legitimately be a decrease in total distance from the scratch, ie a
%positive value isn't just some weird technical issue.... 

%If you're taking the time to read this psuedocode please just feel free to
%mess around with everything :) 


glme_adj = fitglme(modelData, ...
    'DistNorm ~ 1 + TotCells  + scratchArea + genotype*time +  (1|wellPlate) + (1|vidNum)', ...
    'Distribution','normal','Link','identity')




glme_adj_clamp = fitglme(modelData, ...
    'DistNormClamp ~ 1 + TotCells  + scratchArea + genotype*time +  (1|wellPlate) + (1|vidNum)', ...
    'Distribution','normal','Link','identity')




glme_raw = fitglme(modelData, ...
    'DistNorm ~ 1 + TotCells  + scratchArea + genotype*time + (1|wellPlate) + (1|vidNum)', ...
    'Distribution','normal','Link','identity')



normalizeToBaseline=true;
plotStyle='errorBand';% anything other than errorBand will plot all lines
figure
colorMap = struct('KO', [0.1 0.3 0.9], ...      % blue
                  'WT', [19,105,70]/255);       % forest green
alphaVal=.5;
scratchTrackLinePlotterJointAnalysisSimple(structAtPlay, measureAtPlay, ylablr, imFreq, hoursToPlot, ...
                                     genotypes, colorMap, treatmentPrefix, normalizeToBaseline);


% 
% 
% % Raw measures first, normalized next, metadata last
% newOrder = {'DistNorm','startVal','TotCells','scratchArea', ...
%             'Dist','DistNormClamp', ...
%             'time','genotype','vidNum','wellPlate','runID'};
% modelData = modelData(:, newOrder);
% 
% outDir = '/Volumes/PC60/iMicroVids/FinalFigures';
% if ~isfolder(outDir), mkdir(outDir); end
% 
% writetable(modelData, fullfile(outDir, 'RatioInScratchAgg.csv'));   % for Excel/Prism/R
% save(fullfile(outDir, 'RatioInScratchAgg.mat'), 'modelData');       % keeps MATLAB types
% 


%% ok on to single cell 

conversionFactor = 0.65 / 1.13;

trajAtPlay = struct([]);

for rr = 1:size(runDefs,1)
    runID      = runDefs{rr,1};
    geno       = runDefs{rr,2};
    folderPath = runDefs{rr,3};

    if ~isfolder(folderPath), warning('Missing: %s', folderPath); continue, end

    fileList = dir(fullfile(folderPath, '*.mat'));
    fileList = fileList(~startsWith({fileList.name}, '._'));
    fileList = fileList(contains({fileList.name}, 'trajectoryMetrics'));

    for f = 1:numel(fileList)
        loaded = load(fullfile(folderPath, fileList(f).name));
        vs = loaded.velocitySummary;

        tok       = regexp(fileList(f).name, '\((\d+)\)', 'tokens');
        vidNumRaw = str2double(tok{1});
        wellPlate = ceil(vidNumRaw/2);

        vidID  = sprintf('%s_%s_%d', runID, geno, vidNumRaw);
        wellID = sprintf('Well%d_%s', wellPlate, runID);

        for nn = 1:numel(vs)
            rate = vs(nn).Rate_um_per_min;
            if isnumeric(rate) && ~isempty(rate)
                nTP     = numel(rate);
                meanRate = mean(rate*conversionFactor, 'omitnan');
            else
                nTP = 0; meanRate = NaN;
            end

            trajAtPlay(end+1).Genotype  = geno;
            trajAtPlay(end).runID       = runID;
            trajAtPlay(end).vidNum      = vidID;
            trajAtPlay(end).wellPlate   = wellID;
            trajAtPlay(end).NucleusID   = vs(nn).NucleusID;

            trajAtPlay(end).numTimePoints     = nTP;
            trajAtPlay(end).meanRate          = meanRate;
            trajAtPlay(end).NetDisplacement   = vs(nn).NetDisplacement;
            trajAtPlay(end).TotalPath         = vs(nn).TotalPath;
            trajAtPlay(end).displacementRatio = vs(nn).displacementRatio;
            trajAtPlay(end).sharpTurnCount    = vs(nn).sharpTurnCount;
        end
    end
end

fprintf('%d nuclei total\n', numel(trajAtPlay));
disp(groupsummary(struct2table(trajAtPlay(:)), {'runID','Genotype'}))


trajTbl = struct2table(trajAtPlay(:));

cutoff = quantile(trajTbl.TotalPath, 0.8); 

cutoff=158;%    158 is equivalent to 100uM traveled over the duration of the video


trajTbl = trajTbl(   (trajTbl.numTimePoints > 12 & trajTbl.TotalPath> cutoff)   , : )  


trajTbl.Genotype  = reordercats(categorical(trajTbl.Genotype), {'KO','WT'});
trajTbl.runID     = categorical(trajTbl.runID);
trajTbl.wellPlate = categorical(trajTbl.wellPlate);
trajTbl.vidNum    = categorical(trajTbl.vidNum);



colorMap = struct('KO', [0.1 0.3 0.9], ...      % blue
                  'WT', [19,105,70]/255);       % forest green

runIDs   = categories(trajTbl.runID);
runShade = linspace(0, 0.55, numel(runIDs));    % 0 = base color, higher = lighter
genos    = {'KO','WT'};
xPos     = [1 3];

figure; hold on; set(gcf,'Color','w');

for gg = 1:numel(genos)
    geno  = genos{gg};
    color = colorMap.(geno);
    selG  = trajTbl.Genotype == geno;

    % bar = median across all runs for this genotype
    bar(xPos(gg), nanmean(trajTbl.displacementRatio(selG)), ...
        'FaceColor', color, 'FaceAlpha', 0.3, 'BarWidth', 1.2, 'EdgeColor','none');

    % scatter, one shade per run
    for rr = 1:numel(runIDs)
        sel = selG & trajTbl.runID == runIDs{rr};
        n   = sum(sel);
        if n == 0, continue, end

        runColor = color + (1-color)*runShade(rr);
        jitter   = rand(n,1)*0.6 + (xPos(gg) - 0.3);

        scatter(jitter, trajTbl.displacementRatio(sel), 18, ...
                'o', 'MarkerFaceColor', runColor, 'MarkerEdgeColor', 'none', ...
                'MarkerFaceAlpha', 0.6);
    end
end

xlim([0 4])
set(gca, 'XTick', xPos, 'XTickLabel', genos, 'FontSize', 15)
ylabel('Displacement Ratio', 'FontSize', 20)
title('Displacement Ratio by Genotype', 'FontSize', 16)
box on


val2use    = trajTbl.displacementRatio;
genoCell   = cellstr(trajTbl.Genotype);

[p, tbl, stats] = anovan(val2use, {genoCell}, 'varnames', {'Genotype'})

[p, tbl, stats] = anovan(val2use, {cellstr(trajTbl.Genotype), cellstr(trajTbl.wellPlate)},   ...
                         'varnames', {'Genotype','WellPlate'}, 'model', 'linear')





%% single time point Plot
%%%% ---- single time-point bar plot (clamped DistNorm) ----
tPoint = 100;               % minutes; snaps to nearest sampled frame
yVar   = 'DistNorm';   % 'DistNorm' for unclamped
tVals   = unique(modelData.time);
[~, k]  = min(abs(tVals - tPoint));
tActual = tVals(k);
sub     = modelData(modelData.time == tActual & isfinite(modelData.(yVar)), :);

% ---- stats: genotype at this time point, wellPlate as random effect ----
glme_tp = fitglme(sub, sprintf('%s ~ 1 + genotype + (1|wellPlate)', yVar), ...
    'Distribution','normal','Link','identity')
C    = glme_tp.Coefficients;
pTp  = C.pValue(strcmp(C.Name, 'genotype_WT'));

genos  = {'WT','KO'};
runs   = categories(removecats(sub.runID));
nRun   = numel(runs);
shadeF = linspace(-0.35, 0.35, nRun);            % <0 darker, >0 lighter
shade  = @(c,f) (f<0).*(c*(1+f)) + (f>=0).*(c + f*(1-c));
figure; hold on
for gi = 1:numel(genos)
    base = colorMap.(genos{gi});
    yG   = sub.(yVar)(sub.genotype == genos{gi});
    m    = mean(yG);  se = std(yG) / sqrt(numel(yG));
    bar(gi, m, 0.6, 'FaceColor', base, 'FaceAlpha', alphaVal, 'EdgeColor', 'none');
    errorbar(gi, m, se, 'k', 'LineStyle', 'none', 'CapSize', 10, 'LineWidth', 1.2);
    for ri = 1:nRun
        yy = sub.(yVar)(sub.genotype == genos{gi} & sub.runID == runs{ri});
        scatter(gi + (rand(size(yy))-0.5)*0.3, yy, 30, shade(base, shadeF(ri)), ...
                'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
    end
end
yline(0, ':k');
set(gca, 'XTick', 1:numel(genos), 'XTickLabel', genos, 'FontSize', 12, 'Box', 'off');
xlim([0.4 numel(genos)+0.6]);
ylabel('log(value / t_0), negatives set to 0');
title(sprintf('t = %d min, GLME p = %.3g', tActual, pTp));
hold off


