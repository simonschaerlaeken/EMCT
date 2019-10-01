function [ tsvFile ] = emcConvexhull(tsvFile, cfg)
% Computes the convexhull
% 
% syntax
% tsvFile = emcConvexhull(tsvFile, cfg)
% 
% input parameters
% tsvFile: MoCap data structure
% cfg: configuration structure
%     [MANDATORY]
%     -
%     [OPTIONAL]
%     *.featMarker: cell array containing the markernames on which the
%     feature is calculated (default: all)
%     *.frame: int value of the frame at which the body should be displayed
%     *.dispay: boolean deciding if a figure is to be plotted (default: true)
%     if cfg.display: true
%       *.connectionMatrix: [mandatory] list of indexes defining the connection of the
%       body (see emcPlotBody3D)
%     *.dispayUnit: str containing the unit of duration - 's' or 'tf'
% 
%     
% output
% tsvFile: MoCap data structure
% 
% examples
% cfg.featMarker = {'a','b','c'};
% cfg.frame = 40;
% cfg.displayUnit = 's'
% cfg.display = true
% cfg.connectionMatrix = [1 2; 2 3; 1 3]
% tsvFile = emcConvexhull(tsvFile, cfg);
% 
% comments
% feature is saved in tsvFile.processing.convexhull
% 
% see also
% -
% 
% Part of the EMC Toolbox, Copyright 2017,
% University of Geneva, Switzerland
%% CHECKING AREA
if ~isfield(tsvFile, 'info')||~isfield(tsvFile.info, 'filename')
    filename = strsplit(tsvFile.filename, filesep);
    filename = filename{end};
    tsvFile.info.filename = filename;
end
% featMarker - markers on which the feature is calculated
if ~isfield(cfg, 'featMarker')||isempty(cfg.featMarker) 
    disp('cfg.featMarker was not setup. Default: all markers')
    cfg.featMarker = tsvFile.markerName;
end
% display - if true, plot the result
if ~isfield(cfg, 'display')
    disp('cfg.display was not setup. Default: true')
    cfg.display = true;
    errorIfNotField(cfg, 'connectionMatrix')
elseif cfg.display % If true must check for connection matrix
    errorIfNotField(cfg, 'connectionMatrix')
end
% displayUnit - str containing the unit of duration - 's' or 'tf'
if ~isfield(cfg, 'displayUnit')
    disp('cfg.displayUnit was not setup. Default: "tf"')
    cfg.displayUnit = 'tf';
end
% frame - frame at which the body will be plotted
if ~isfield(cfg, 'frame')||isempty(cfg.frame) % Must choose random frame
    % Generate random point in time
    a = 1;
    b = size(tsvFile.data,1);
    cfg.frame = int64(round((b-a).*rand(1,1) + a));    
end


%% COMPUTATION AREA
disp('[PROCESSING] Convexhull')
% Select markers
featMarkerIndex = findIndexList(tsvFile.markerName, cfg.featMarker);
tsvFileMarker = mcgetmarker(tsvFile, featMarkerIndex);
convexhull = zeros(size(tsvFileMarker.data,1),1);
% ft_c_hull = zeros(size(tsv_file.data,1),1);
% Go through the each timeframe
for i = 1:size(tsvFileMarker.data,1)
    % Create list of points
    listPoints = (reshape(tsvFileMarker.data(i,:),3,length(tsvFileMarker.data(i,:))/3))';
    [~,v] = convhull(listPoints(:,1),listPoints(:,2),listPoints(:,3));
    convexhull(i) = v;
end
if cfg.display
    figure('Name',[tsvFile.info.filename, ' - Convexhull'],'NumberTitle','off')
    if strcmp(cfg.displayUnit, 's')
        time = 0:(1/tsvFile.freq):((tsvFile.nFrames-1)/tsvFile.freq);
        plot(time, convexhull)
    else
        plot(convexhull)
    end
end

if cfg.display
    listPoints = (reshape(tsvFileMarker.data(cfg.frame,:),3,length(tsvFileMarker.data(cfg.frame,:))/3))';
    [k,~] = convhull(listPoints(:,1),listPoints(:,2),listPoints(:,3));
    % Plot
    figure('Name',[tsvFile.info.filename, ' - Convexhull - TF=', num2str(cfg.frame)],'NumberTitle','off')
    h1 = subplot(1,2,1);
    if strcmp(cfg.displayUnit, 's')
        time = 0:(1/tsvFile.freq):((tsvFile.nFrames-1)/tsvFile.freq);
        plot(time, convexhull)
    else
        plot(convexhull)
    end
    title(h1,'Volume of Convex Hull over time')
    h2 = subplot(1,2,2);
    s = trisurf(k,listPoints(:,1),listPoints(:,2),listPoints(:,3),'Facecolor','cyan');
    alpha(s, .1);
    title(h2, 'Convex Hull at a specific point in time')
    hold on;
    % Superimpose the skeleton on the convex hull
    axis equal;
    set(gcf, 'renderer', 'opengl');
    set(gca, 'CameraPosition', [400 -200 800]);
    % Create all the points
    for i = 1:tsvFile.nMarkers
        markerIn = i*3-2;
        markerOut = markerIn+2;
        points(i,:) = tsvFile.data(cfg.frame,markerIn:markerOut);
    end
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
        drawEdge(edge, 'color', 'r', 'linewidth', 2);
    end
end

% Return the final value
tsvFile.processing.convexhull = convexhull;
end


