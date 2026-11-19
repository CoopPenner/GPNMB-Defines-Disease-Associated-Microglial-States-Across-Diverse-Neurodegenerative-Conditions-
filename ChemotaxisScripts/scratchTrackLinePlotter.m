function [modelData] = scratchTrackLinePlotter(structAtPlay,measureAtPlay, ylablr, imFreq,hoursToPlot, genotypes,colorMap,alphaVal,treatmentPrefix,normalizeToBaseline,plotStyle)

figure;

%  containers
traceCatchTot = [];
trace2PlotCatchTot=[];
genoCatchTot = [];
timeCatchTot = [];
vidNumCatchTot = {};
countCatchTot = [];
startValCatchTot=[];
scratchAreaCatchTot = [];
wellPlateCatchTot=[];
traceMat=[];
maxLength = hoursToPlot*(60/imFreq); % 6 hours for first set captured up to 14 hours
maxLength=maxLength-1; %there is sometimes a slight mismatch so I just cut the last timepoint
timeVec =   [0:imFreq:imFreq* (maxLength-1)] ;% Time in minutes

for gg = 1:numel(genotypes)
    geno = genotypes{gg};
    color = colorMap.(geno);

    % Get all entries for this genotype
    genoEntries = structAtPlay(strcmp({structAtPlay.Genotype}, geno));
    
    % Extract vectors
    if strcmp(measureAtPlay,'globalDist')
    traces = {genoEntries.AvgDistToScratch};
    elseif strcmp(measureAtPlay, 'scratchCount')
        traces={genoEntries.NearScratchCount};
    end
    totCounts = {genoEntries.TotObjects};
    vidNums = [genoEntries.VideoNum];
    wellPlate=[genoEntries.wellPlate];
    scratchAreas = [genoEntries.ScratchArea];

    % subplot(1, 2, gg);
    hold on;

    for jj = 1:numel(traces)


colorAtPlay=color;

if wellPlate(jj)==1
    colorAtPlay(1)=.25;
elseif wellPlate(jj)==2
    colorAtPlay(1)=.5;
elseif wellPlate(jj)==3
    colorAtPlay(1)=.75;
elseif wellPlate(jj)==4
    colorAtPlay(1)=1;
end

        
        trace = traces{jj}(1:maxLength);
        countz=totCounts{jj}(1:maxLength);


        if strcmp(measureAtPlay, 'scratchCount')
        trace = traces{jj}(1:maxLength);
        countz=totCounts{jj}(1:maxLength);

        trace=trace'./countz';



        end



        if normalizeToBaseline
    trace2Plot=trace./ trace(1);
    trace2PlotCatchTot=[trace2PlotCatchTot,trace2Plot'];
        else
            trace2Plot=trace;
            trace2PlotCatchTot=[trace2PlotCatchTot,trace2Plot'];


end



if ~strcmp(plotStyle, 'errorBand')
        % Plot individual trace
        plot(timeVec, trace2Plot, 'Color', [colorAtPlay alphaVal], 'LineWidth', 1);

end
 

        
        % Store values for GLME
        traceCatchTot = [traceCatchTot, trace'];
        genoCatchTot = [genoCatchTot, repmat({geno}, 1, maxLength)];
        timeCatchTot = [timeCatchTot, timeVec ];
        vidNumCatchTot = [vidNumCatchTot, repmat({sprintf('%s_%d', geno, vidNums(jj))}, 1, maxLength)];
        wellPlateCatchTot=[wellPlateCatchTot,repmat({sprintf('%s_%d', geno, wellPlate(jj))}, 1, maxLength) ];
        countCatchTot = [countCatchTot,countz' ];
        scratchAreaCatchTot = [scratchAreaCatchTot, (scratchAreas(jj)* ones(1, maxLength))];
        startValCatchTot=[startValCatchTot, (trace(1)* ones(1, (maxLength)))];


    end

    % Plot mean trace


   % allTraces = cellfun(@(x) x(1:maxLength), traces, 'UniformOutput', false);
    
    %traceMat = cell2mat(cellfun(@(x) padarray(x(:)', [0, maxLength - length(x)], NaN, 'post'), allTraces, 'UniformOutput', false)');
    
    % reshape into rows = traces, columns = timepoints
        traceMat = reshape(trace2PlotCatchTot, maxLength, [])';


        traceMat(isinf(traceMat))=nan;
        
        meanTrace = nanmean(traceMat, 1);

       plot(timeVec, meanTrace, '-', 'Color', color, 'LineWidth', 2.5);

 nPerT     = sum(~isnan(traceMat), 1);
 semTrace  = nanstd(traceMat, 0, 1) ./ sqrt(nPerT);      % for SD instead: nanstd(traceMat,0,1)

    % shaded error band (errorband mode only)
    if strcmp(plotStyle, 'errorBand')
        xFill = [timeVec, fliplr(timeVec)];
        yFill = [meanTrace + semTrace, fliplr(meanTrace - semTrace)];
        fill(xFill, yFill, color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    end



    % Format subplot
    xlabel('Time (Minutes)','FontSize',20);
    ylabel(ylablr,'FontSize', 20);
   % title(sprintf('%s Microglia Treated with %s', geno, treatmentPrefix),'FontSize',30);
    box on;
    xlim([0, maxLength*imFreq]);






 %        subplot(2,2,3:4)
 %        hold on
 %    plot(timeVec, meanTrace, '-', 'Color', color, 'LineWidth', 2.5);
 % xlabel('Time (Minutes)');
 %    ylabel('Average Distance to Scratch');
 %    title(sprintf('%s Microglia Treated with %s', geno, treatmentPrefix));
 %    box on;
 %    xlim([0, maxLength*frameRate]);


%statistical Testing using glme


% make model  table
modelData = table( ...
    traceCatchTot(:), ...
    startValCatchTot(:),...
    vidNumCatchTot(:), ...
    countCatchTot(:), ...
    scratchAreaCatchTot(:), ...
    genoCatchTot(:), ...
    timeCatchTot(:), ...
    wellPlateCatchTot(:),...
    'VariableNames', {'Dist', 'startVal', 'vidNum', 'TotCells', 'scratchArea', 'genotype', 'time', 'wellPlate'} ...
);


modelData.vidNum = categorical(modelData.vidNum);  % Keep only if each entry is a scalar
modelData.genotype = categorical(modelData.genotype);
modelData.wellPlate= categorical(modelData.wellPlate);


end
