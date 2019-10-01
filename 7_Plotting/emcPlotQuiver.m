function [] = emcPlotQuiver(tsvFile, cfg)
% Plots markers, their line of movement in
% space as well as arrows for Speed and acceleration using the function
% Quiver from Matlab
% 
% syntax
% emcPlotQuiver(tsvFile, cfg);
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     -
%     [OPTIONAL]
%     *.step: int value defining the step between two arrows (default: 20)
%     *.plotMarkerList: cell array with markers to be plotted (if not defined all markers will be used)
%     *.plotBodyFlag: boolean flag defining if the body structure is
%     printed or not
%     if cfg.plotBodyFlag: true
%       *.connectionMatrix: list of indexes defining the connection of the
%       body (see emcPlotBody3D)
%     
% output
% -
% 
% examples
% cfg.step = 40;
% cfg.plotMarkerList = {'a','b','c'};
% cfg.plotBodyFlag = true;
% cfg.connectionMatrix = [1 2; 2 3; 1 3]
% emcPlotQuiver(tsvFile, cfg);
% 
% comments
% -
% 
% see also
% emcPlotBody3D
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
% step - step at which the points will be displayed
if ~isfield(cfg, 'step')||isempty(cfg.step)
    disp('INFO: Since no cfg.step was set up, arbirarily defined to 20.')
    cfg.step = 20;
end
% plotMarkerList - markers that will be displayed
if ~isfield(cfg, 'plotMarkerList')||isempty(cfg.plotMarkerList)
    disp('INFO: Since no cfg.plotMarkerList was set up, arbirarily defined to all markers.')
    cfg.plotMarkerList = tsvFile.markerName;
end
% plotBodyFlag - Flag, if true, will plot the structure of the body onto
% the graph
if ~isfield(cfg, 'plotBodyFlag')
    cfg.plotBodyFlag = true;
end
if cfg.plotBodyFlag
    errorIfNotField(cfg, 'connectionMatrix')
    cfg.frame = 1;
end


%% COMPUTATION AREA
disp(['[PLOT] ', tsvFile.info.filename ' - step: ' num2str(cfg.step)]);
% Extract data using just the markers we want to compare
listMarkersSelectedIdx = findIndexList(tsvFile.markerName,cfg.plotMarkerList);
tsvFileMarker = mcgetmarker(tsvFile,listMarkersSelectedIdx);
% Derivation
warning('off','all')
tsvFileSpeed = mctimeder(tsvFileMarker, 1);
tsvFileAcc = mctimeder(tsvFileMarker, 2);
warning('on','all')
figure('Name',[tsvFile.info.filename ' - step: ' num2str(cfg.step)],'NumberTitle','off')
% Plot
for markerIdx = 1:length(cfg.plotMarkerList)
    hold on
    indexStart = (((markerIdx)-1)*3)+1;
    plot3(tsvFileMarker.data(:,indexStart),... % Trajectory
        tsvFileMarker.data(:,indexStart+1),...
        tsvFileMarker.data(:,indexStart+2),'k')
    quiver3(tsvFileMarker.data(1:cfg.step:size(tsvFileMarker.data,1),indexStart),... % Speed
        tsvFileMarker.data(1:cfg.step:size(tsvFileMarker.data,1),indexStart+1),...
        tsvFileMarker.data(1:cfg.step:size(tsvFileMarker.data,1),indexStart+2),...
        tsvFileSpeed.data(1:cfg.step:size(tsvFileMarker.data,1),indexStart),...
        tsvFileSpeed.data(1:cfg.step:size(tsvFileMarker.data,1),indexStart+1),...
        tsvFileSpeed.data(1:cfg.step:size(tsvFileMarker.data,1),indexStart+2), 'r')
    quiver3(tsvFileMarker.data(1:cfg.step:size(tsvFileMarker.data,1),indexStart),... % Acceleration
        tsvFileMarker.data(1:cfg.step:size(tsvFileMarker.data,1),indexStart+1),...
        tsvFileMarker.data(1:cfg.step:size(tsvFileMarker.data,1),indexStart+2),...
        tsvFileAcc.data(1:cfg.step:size(tsvFileMarker.data,1),indexStart),...
        tsvFileAcc.data(1:cfg.step:size(tsvFileMarker.data,1),indexStart+1),...
        tsvFileAcc.data(1:cfg.step:size(tsvFileMarker.data,1),indexStart+2), 'g')
    legend({'Trajectory','Speed','Acceleration'})
end
for markerIdx = 1:length(cfg.plotMarkerList)
    indexStart = (((markerIdx)-1)*3)+1;
    drawPoint3d(tsvFileMarker.data(1,indexStart:indexStart+2), 'marker', 'o', 'markersize', 8, 'linewidth', 2, 'markerFaceColor', 'w');
    drawPoint3d(tsvFileMarker.data(end,indexStart:indexStart+2), 'marker', 'x', 'markersize', 8, 'linewidth', 2, 'markerFaceColor', 'w');
end

% Print structure
if isfield(cfg, 'plotBodyFlag') && cfg.plotBodyFlag
    % Create all the points
    points = zeros(tsvFile.nMarkers,3);
    for i = 1:tsvFile.nMarkers,
        markerIn = i*3-2;
        markerOut = markerIn+2;
        points(i,:) = tsvFile.data(cfg.frame,markerIn:markerOut);
    end
    % prepare a figure for drawing
    axis equal;
    set(gcf, 'renderer', 'opengl');
    set(gca, 'CameraPosition', [400 -200 800]);
    % Draw all the edges
    for i = 1:size(cfg.connectionMatrix,1),
        edge = zeros(1,6);
        for j = 1:size(cfg.connectionMatrix,2),
            markerInEdge = j*3-2;
            markerOutEdge = markerInEdge+2;
            edge(markerInEdge:markerOutEdge) = points(cfg.connectionMatrix(i,j),:);
        end
        drawEdge(edge, 'color', 'k', 'linewidth', 1);
    end
end

hold off
view([45,45,45]);
end
