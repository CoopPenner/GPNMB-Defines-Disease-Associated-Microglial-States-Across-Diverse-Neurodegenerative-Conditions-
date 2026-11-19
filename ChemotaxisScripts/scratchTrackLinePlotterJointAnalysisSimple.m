function scratchTrackLinePlotterJointAnalysisSimple(structAtPlay, measureAtPlay, ylablr, imFreq, hoursToPlot, ...
                                 genotypes, colorMap, treatmentPrefix, normalizeToBaseline)

figure; hold on; %set(gcf,'Color','w');

maxLength = hoursToPlot*(60/imFreq) - 1;
timeVec   = 0:imFreq:imFreq*(maxLength-1);

legH = []; legL = {};

for gg = 1:numel(genotypes)
    geno  = genotypes{gg};
    color = colorMap.(geno);

    sel = structAtPlay(strcmp({structAtPlay.Genotype}, geno));
    if isempty(sel), continue, end

    switch measureAtPlay
        case 'globalDist',   traces = {sel.AvgDistToScratch};
        case 'scratchCount', traces = {sel.NearScratchCount};
    end
    totCounts = {sel.TotObjects};

    genoTraceMat = [];
    for jj = 1:numel(traces)
        trace  = traces{jj}(1:maxLength);
        countz = totCounts{jj}(1:maxLength);

        if strcmp(measureAtPlay,'scratchCount')
            trace = trace(:)./countz(:);
        end
        trace = trace(:)';

        if normalizeToBaseline
            trace = trace./trace(1);
        end
        trace(isinf(trace)) = NaN;

        genoTraceMat(end+1,:) = trace;
    end

    meanTrace = nanmean(genoTraceMat, 1);
    nPerT     = sum(~isnan(genoTraceMat), 1);
    semTrace  = nanstd(genoTraceMat, 0, 1) ./ sqrt(nPerT);



fill([timeVec fliplr(timeVec)], ...
     [meanTrace+semTrace fliplr(meanTrace-semTrace)], ...
     color, 'FaceAlpha', 0.2, 'EdgeColor','none');set(gcf,'Renderer','painters');

    h = plot(timeVec, meanTrace, '-', 'Color', color, 'LineWidth', 3);
    legH(end+1) = h; legL{end+1} = sprintf('%s (n=%d)', geno, size(genoTraceMat,1));
end

xlabel('Time (Minutes)','FontSize',20);
ylabel(ylablr,'FontSize',20);
xlim([0, maxLength*imFreq]);
%legend(legH, legL, 'Location','best','Box','off');
box on;
end