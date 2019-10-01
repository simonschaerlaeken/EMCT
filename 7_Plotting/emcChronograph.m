function [] = emcChronograph(tsvFile, cfg)
% Plots a sequence of the full movement frame after frame in one image
% (chronophotographic) and a boxed representation of the magnitude of the
% defined feature
% 
% syntax
% emcVisualSequence(tsvFile, cfg);
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     *.connectionMatrix: list of indexes defining the connection of the
%       body (see emcPlotBody3D)
%     *.visSeqFeature: str feature name to be displayed
%     [OPTIONAL]
%     *.step: int number of frames to skip in between two timestamps (default: 50)
%     *.width: str width of a single frame within the chronophotographic (default: 40)
%     
% output
% -
% 
% examples
% cfg.connectionMatrix = [1 2; 2 3; 1 3]
% cfg.visSeqFeature = 'speed';
% cfg.step = 30;
% cfg.width = 50;
% emcVisualSequence(tsvFile, cfg);
% 
% 
% comments
% -
% 
% see also
% -
% 
% Part of the EMC Toolbox, Copyright 2017,
% University of Geneva, Switzerland
%% CHECKING AREA
if ~isfield(tsvFile, 'info')||~isfield(tsvFile.info, 'filename'),
    filename = strsplit(tsvFile.filename, filesep);
    filename = filename{end};
    tsvFile.info.filename = filename;
end
% Error if missing
% step - define the number of frames to skip in between two timestamps.
if ~isfield(cfg, 'step')
    cfg.step = 50;
end

% width - define the width of a single frame within the chronophotographic.
if ~isfield(cfg, 'width')
    cfg.width = 40;
end

% visSeqFeature - define the feature to print.
if ~isfield(cfg, 'visSeqFeature')
    cfg.visSeqFeature = '';
else
    if ~isfield(tsvFile, 'processing')
        error('processing field not found in TSVFILE, check if the processing step was computed')
    else
        if ~isfield(tsvFile.processing, cfg.visSeqFeature)
            error([cfg.visSeqFeature, ' is not part of the computed feature for this TSV FILE'])
        end
    end
 
end


%% COMPUTATION AREA
% Info for window
nbSeq = size(tsvFile.data,1)/cfg.step;
heightTotal = 400;
widthTotal = nbSeq*cfg.width/3;
% Color map
colorMap = hot(8); % spring, winter, summer,...
% Draw fig
hVisSeq = figure('Name',[tsvFile.info.filename, ' - Visual Sequence - ', cfg.visSeqFeature],'NumberTitle','off');
set(hVisSeq, 'Position', [0 100 widthTotal heightTotal])
hGraph = gca;
if ~isempty(cfg.visSeqFeature)
    % Feat min max
    % -- Average accros markers
    feature = tsvFile.processing.(cfg.visSeqFeature);
    featureAvg = mean(feature,2);
    % -- Min Max
    featureMin = min(featureAvg);
    featureMax = max(featureAvg);
end
% Init
lowestPoint = findLowestPoint(tsvFile); % Start high and then gets replaced by min values
shiftIdx = 0;
colorShift = 0.5;
colorIncrement = (1-colorShift)/nbSeq;
for frame = 1:cfg.step:size(tsvFile.data,1), 
    % Create all the points
    points = zeros(tsvFile.nMarkers,2);
    for i = 1:tsvFile.nMarkers,
        markerIn = i*3-2;
        markerOut = markerIn+2;
        points(i,:) = [tsvFile.data(frame,markerIn),tsvFile.data(frame,markerOut)];
    end
    % Shift them depending on the frame according to X
    points(:,1) = points(:,1) + shiftIdx*cfg.width;
    % Color
    colorBody = [1,1,1]*colorShift;
    % Print body
    plotBody2DColor(hGraph, cfg, points, colorBody) 
    % Print Rectangle
    if ~isempty(cfg.visSeqFeature)
        featurePower = featurePowerCompute(cfg, featureAvg, frame, featureMin, featureMax);
        colorMapIdx = ceil(featurePower*size(colorMap,1)); % Ceil instead of round to have at least idx = 1
        colorRect = colorMap(colorMapIdx,:);
        drawrectangle(hGraph, cfg, shiftIdx, colorRect, lowestPoint)
    end
    % Next
    shiftIdx = shiftIdx + 1;
    colorShift = colorShift - colorIncrement;
end
% AXIS
% Remove Y ticks
hGraph.YTick = [];
% Change X ticks
hGraph.XTick = 0:cfg.width:nbSeq*cfg.width;
hGraph.XTickLabel  = 1:cfg.step:nbSeq*cfg.step;
% Color axis scaling
if ~isempty(cfg.visSeqFeature)
    caxis([featureMin featureMax])
end
end



function [] = plotBody2DColor(hGraph, cfg, points, color)   
% Function to plot in 3D the representation of the point and lines
% ----------------------------------------------------------------
% Plot
hold(hGraph, 'on');
set(gcf, 'renderer', 'opengl');
% Draw all points
drawPoint(hGraph,points, 'color', color);
% Draw all the edges
for i = 1:size(cfg.connectionMatrix,1),
    edge = zeros(1,4);
    for j = 1:size(cfg.connectionMatrix,2),
        markerInEdge = j*2-1;
        markerOutEdge = markerInEdge+1;
        edge(markerInEdge:markerOutEdge) = points(cfg.connectionMatrix(i,j),:);
    end
    drawEdge(hGraph,edge, 'color', color, 'linewidth', 1);
end
hold(hGraph, 'off');
end

function [featurePower] = featurePowerCompute(cfg, featureAvg, frame, featureMin, featureMax)
% Function which returns the average value of a feature for a specific
% point in the sequence (average around it by the framskip, and over all markers)

% Select Chunk
startIdx = round(frame-cfg.step/2);
if startIdx < 1,
    startIdx = 1;
end
endIdx = round(frame+cfg.step/2);
if endIdx > size(featureAvg,1),
    endIdx = size(featureAvg,1);
end
featureChunk = featureAvg(startIdx:endIdx,:);

% Average over time
featureAvgAvg = mean(featureChunk,1);

% Compute feature power
featurePower = (featureAvgAvg-featureMin)/(featureMax-featureMin);
end

function [] = drawrectangle(hGraph, cfg, shiftIdx, color, lowestPoint)
% Function to draw a rectangle with a specific width equal to the frame
% skipped, and a specific color corresponding to the power of the feature
% Plot
hold(hGraph, 'on');
width = cfg.width;
height = 200;
xpos = (shiftIdx*cfg.width)-(cfg.width/2); % centered around the position of the body
ypos = lowestPoint-height-50; % -50 is to decrease
rectangle('Position',[xpos,ypos,width,height],'FaceColor', color,'Parent', hGraph);
hold(hGraph, 'off');
end

function lowestPoint = findLowestPoint(tsv_file)
% Function which find the lowest point in a csv (on the Z dimension)
data = tsv_file.data(3:3:end);
lowestPoint = min(min(data));
end

