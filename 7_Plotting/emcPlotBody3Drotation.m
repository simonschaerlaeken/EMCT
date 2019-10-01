function [] = emcPlotBody3D(tsvFile, cfg)
% Plots the structure of the body in 3D
% 
% syntax
% emcPlotBody3D(tsvFile, cfg);
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     *.connectionMatrix: list of indexes defining the connection of the
%       body 
%     [OPTIONAL]
%     *.frame: int value defining the frame at which the marker should be
%     printed (default: random)
%     *.lineVector: cell array containing pairs of markers defining red
%     lines
%     
% output
% -
% 
% examples
% cfg.connectionMatrix = [1 2; 2 3; 1 3]
% cfg.frame = 40;
% cfg.lineVector = {{'a','b'},{'b','c'}};
% emcPlotBody3D(tsvFile, cfg);
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
% filename - for title
if ~isfield(tsvFile, 'info')||~isfield(tsvFile.info, 'filename')
    filename = strsplit(tsvFile.filename, filesep);
    filename = filename{end};
    tsvFile.info.filename = filename;
end
% frame - frame at which the body will be plotted
if ~isfield(cfg, 'frame')||isempty(cfg.frame)||cfg.frame>size(tsvFile.data,1)
    % Must choose random frame
    % Generate random point in time
    a = 0;
    b = size(tsvFile.data,1);
    cfg.frame = int64(round((b-a).*rand(1,1) + a));  
end
% connectionMatrix - table containing all the connexions (sticks) between
% markers. Necessary to run
errorIfNotField(cfg, 'connectionMatrix')

errorIfNotField(tsvFile, 'rotationdata')

% lineVector - cell list containing marker to make lines in red
if ~isfield(cfg, 'lineVector')
    cfg.lineVector = {};
elseif numel(cfg.lineVector) == 2 && ischar(cfg.lineVector{1}) % Only one line
    cfg.lineVector = {cfg.lineVector};
else
    for lineIdx = 1:numel(cfg.lineVector)
        if numel(cfg.lineVector{lineIdx}) ~= 2
            error('cfg.lineVector must be of this format: {{"a","b"},{"v","d"}}')
        end
    end
end
%% COMPUTATION AREA
disp(['[PLOT] ', tsvFile.info.filename ' - frame: ' num2str(cfg.frame)]);
% Create all the points
points = zeros(tsvFile.nMarkers,3);
for i = 1:tsvFile.nMarkers
    markerIn = i*3-2;
    markerOut = markerIn+2;
    points(i,:) = tsvFile.data(cfg.frame,markerIn:markerOut);
end
% prepare a figure for drawing
figure('Name',[tsvFile.info.filename ' -- TF = ' num2str(cfg.frame)],'NumberTitle','off'); clf; hold on;
axis equal;
set(gcf, 'renderer', 'opengl');
set(gca, 'CameraPosition', [400 -200 800]);
% Draw all points
drawPoint3d(points);
% Draw all the edges
for i = 1:size(cfg.connectionMatrix,1)
    edge = zeros(1,6);
    for j = 1:size(cfg.connectionMatrix,2)
        markerInEdge = j*3-2;
        markerOutEdge = markerInEdge+2;
        edge(markerInEdge:markerOutEdge) = points(cfg.connectionMatrix(i,j),:);
    end
    drawEdge(edge, 'color', 'k', 'linewidth', 1);
end

if ~isempty(cfg.lineVector)
    for i = 1:numel(cfg.lineVector)
        lineMarkers = cfg.lineVector{i};
        markersIdx = findIndexList(tsvFile.markerName, lineMarkers); % Find the 2 markers creating the vectors
        line = createLine3d(points(markersIdx(1),:), points(markersIdx(2),:));
        drawLine3d(line, 'color', 'r','linewidth', 2);
    end
end
% Draw Rotation
for markerIdx = 1:numel(tsvFile.rotationMarkerName)
    posIdx = ((markerIdx-1)*3)+1;
    rotIdx = ((markerIdx-1)*tsvFile.rotationDim)+1;
    quiver3( tsvFile.data(cfg.frame,posIdx),tsvFile.data(cfg.frame,posIdx+1),tsvFile.data(cfg.frame,posIdx+2),...
             tsvFile.rotationdata(cfg.frame,rotIdx+0),tsvFile.rotationdata(cfg.frame,rotIdx+1),tsvFile.rotationdata(cfg.frame,rotIdx+2),'r');
    quiver3( tsvFile.data(cfg.frame,posIdx),tsvFile.data(cfg.frame,posIdx+1),tsvFile.data(cfg.frame,posIdx+2),...
             tsvFile.rotationdata(cfg.frame,rotIdx+3),tsvFile.rotationdata(cfg.frame,rotIdx+4),tsvFile.rotationdata(cfg.frame,rotIdx+5),'g');
    quiver3( tsvFile.data(cfg.frame,posIdx),tsvFile.data(cfg.frame,posIdx+1),tsvFile.data(cfg.frame,posIdx+2),...
             tsvFile.rotationdata(cfg.frame,rotIdx+6),tsvFile.rotationdata(cfg.frame,rotIdx+7),tsvFile.rotationdata(cfg.frame,rotIdx+8),'b');
end
hold off;

end

