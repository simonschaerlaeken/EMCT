function [] = emcPlotBody2D(tsvFile, cfg)
%EMCPLOTBODY2D Function to plot the structure of the body in 2D
%   IN: tsvFile - structure defined as TSV
%       cfg     - configuration structure containing the options

%% CHECKING AREA
% filename - for title
if ~isfield(tsvFile, 'info')||~isfield(tsvFile.info, 'filename'),
    filename = strsplit(tsvFile.filename, filesep);
    filename = filename{end};
    tsvFile.info.filename = filename;
end
% frame - frame at which the body will be plotted
if ~isfield(cfg, 'frame')||isempty(cfg.frame),
    % Must choose random frame
    % Generate random point in time
    a = 1;
    b = size(tsvFile.data,1);
    cfg.frame = int64(round((b-a).*rand(1,1) + a));  
end
% connectionMatrix - table containing all the connexions (sticks) between
% markers. Necessary to run
errorIfNotField(cfg, 'connectionMatrix')

% lineVector - cell list containing marker to make lines in red
if ~isfield(cfg, 'lineVector')
    cfg.lineVector = {};
elseif numel(cfg.lineVector) == 2 && ischar(cfg.lineVector{1}) % Only one line
    cfg.lineVector = {cfg.lineVector};
else
    for lineIdx = 1:numel(cfg.lineVector),
        if numel(cfg.lineVector{lineIdx}) ~= 2
            error('cfg.lineVector must be of this format: {{"a","b"},{"v","d"}}')
        end
    end
end

%% COMPUTATION AREA
disp(['[PLOT] ', tsvFile.info.filename ' - frame: ' num2str(cfg.frame)]);
% Create all the points
points = zeros(tsvFile.nMarkers,2);
for i = 1:tsvFile.nMarkers,
    markerIn = i*2-1;
    markerOut = markerIn+1;
    points(i,:) = tsvFile.data(frame,markerIn:markerOut);
end
% prepare a figure for drawing
figure('Name',[tsvFile.info.filename ' -- TF = ' num2str(cfg.frame)],'NumberTitle','off'); clf; hold on;
axis equal;
set(gcf, 'renderer', 'opengl');
hold on
% Draw all points
drawPoint(points);
% Draw all the edges
for i = 1:size(cfg.connectionMatrix,1),
    edge = zeros(1,4);
    for j = 1:size(cfg.connectionMatrix,2),
        markerInEdge = j*2-1;
        markerOutEdge = markerInEdge+1;
        edge(markerInEdge:markerOutEdge) = points(cfg.connectionMatrix(i,j),:);
    end
    drawEdge(edge, 'color', 'k', 'linewidth', 1);
end

if ~isempty(cfg.lineVector)
    for i = 1:numel(cfg.lineVector),
        lineMarkers = cfg.lineVector{i};
        markersIdx = findIndexList(tsvFile.markerName, lineMarkers); % Find the 2 markers creating the vectors
        line = createLine(points(markersIdx(1),:), points(markersIdx(2),:));
        drawLine(line, 'color', 'r','linewidth', 2);
    end
end

hold off
end

